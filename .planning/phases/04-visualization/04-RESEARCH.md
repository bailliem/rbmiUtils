# Phase 4: Visualization - Research

**Researched:** 2026-02-08
**Domain:** ggplot2 forest plots for clinical trial treatment effects
**Confidence:** HIGH

## Summary

This phase implements `plot_forest()`, a function that takes an rbmi pool object and produces a publication-quality forest plot with an aligned table panel. The function builds on the existing `tidy_pool_obj()` infrastructure (Phase 1) to extract treatment effects and/or LSM estimates, then renders them as a horizontal forest plot using ggplot2 with a patchwork-composed table panel.

The standard approach in clinical R packages for forest plots with aligned tables is to create multiple ggplot2 objects (a table panel using `geom_text()` + `theme_void()`, and the forest plot itself using `geom_linerange()` + `geom_point()`) and combine them via patchwork with controlled `widths`. This pattern is used by KHstats, tern/NEST (Roche), and multiple pharma visualization workflows. The Okabe-Ito palette (built into base R since 4.0) is the established colorblind-friendly standard for clinical visualizations.

**Primary recommendation:** Use a patchwork-based three-panel approach (table | forest plot | p-values) with ggplot2 geom_linerange + geom_point, returning the composed patchwork object (which inherits from ggplot and supports `& theme()` for user customization).

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| ggplot2 | 4.0.1+ | Plot engine | Already in Suggests; the only viable ggplot-based plot system |
| patchwork | 1.3.2+ | Multi-panel composition | De facto standard for combining ggplots; supports `& theme()` for user override |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| grDevices (base R) | built-in | `palette.colors("Okabe-Ito")` | Color palette; no external dependency needed |
| scales | 1.4.0+ | Axis formatting (already a ggplot2 dependency) | If custom axis label formatting needed |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| patchwork | cowplot::plot_grid() | cowplot is installed but patchwork has better width control and native ggplot inheritance |
| patchwork | gridExtra::grid.arrange() | gridExtra is less ergonomic; no `& theme()` support |
| patchwork | Pure ggplot2 clip="off" with margin | Fragile; text width varies by content; no reliable alignment |
| ggplot2 | forestploter | forestploter is not ggplot2-based; users cannot customize with `+ theme()` |
| ggplot2 | forestly | Interactive/HTML only; not suitable for CSR/PDF output |

**Dependency strategy:**
- ggplot2: Already in `Suggests` -- keep there; guard with `requireNamespace()` at function entry (same pattern as `efficacy_table()` with gt)
- patchwork: Add to `Suggests` -- guard with `requireNamespace()` at function entry
- No new `Imports` dependencies needed

**Installation (for users):**
```r
install.packages(c("ggplot2", "patchwork"))
```

## Architecture Patterns

### Recommended Project Structure
```
R/
  plot_forest.R          # Main plot_forest() function + internal helpers
  plot_utils.R           # Shared plot utilities (theme, palette, etc.)
tests/testthat/
  test-plot_forest.R     # Tests for plot_forest()
```

Alternatively, `plot_utils.R` content could live inside `plot_forest.R` if the helpers are minimal. Given that this is the only visualization function (per phase scope), keeping it in a single file is cleaner.

### Pattern 1: Three-Panel Patchwork Composition
**What:** Build the forest plot as three aligned ggplot2 objects combined via patchwork.
**When to use:** Always -- this is the recommended approach for forest plots with table panels.
**Structure:**
```
[Left: Table Panel]  |  [Middle: Forest Plot]  |  [Right: P-value Panel]
  Visit labels          Point + CI whiskers       P-value text
  Estimate (95% CI)     Reference line at 0
```

**Example (verified on installed ggplot2 4.0.1 + patchwork 1.3.2):**
```r
# Left panel: visit labels + estimate text
p_left <- ggplot(plot_data, aes(y = visit)) +
  geom_text(aes(x = 0, label = est_ci_label), hjust = 0.5, size = 3) +
  scale_y_discrete(limits = rev) +
  theme_void() +
  theme(axis.text.y = element_text(hjust = 1, size = 9))

# Middle panel: forest plot
p_mid <- ggplot(plot_data, aes(x = est, y = visit, xmin = lci, xmax = uci)) +
  geom_linerange(linewidth = 0.6) +           # CI whiskers (no caps)
  geom_point(size = 3, shape = 16) +           # Filled circles
  geom_vline(xintercept = ref_value, linewidth = 0.5) +  # Solid reference line
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )

# Right panel: p-values
p_right <- ggplot(plot_data, aes(y = visit)) +
  geom_text(aes(x = 0, label = pval_label), hjust = 0.5, size = 3) +
  scale_y_discrete(limits = rev) +
  theme_void()

# Combine
combined <- p_left + p_mid + p_right +
  patchwork::plot_layout(widths = c(3, 4, 1.5))
```

### Pattern 2: Data Preparation Pipeline (Reuse tidy_pool_obj)
**What:** Use existing `tidy_pool_obj()` to parse the pool object, then filter/transform for plotting.
**When to use:** Always -- consistent with `efficacy_table()` pattern.
**Example:**
```r
plot_forest <- function(pool_obj, ...) {
  # Step 1: Tidy the pool object (reuse Phase 1 infrastructure)
  tidy_df <- tidy_pool_obj(pool_obj)

  # Step 2: Filter to requested display mode
  if (display == "trt") {
    plot_data <- tidy_df[tidy_df$parameter_type == "trt", ]
  } else {
    plot_data <- tidy_df[tidy_df$parameter_type == "lsm", ]
  }

  # Step 3: Clean visit labels (same as efficacy_table pattern)
  visit_clean <- gsub("_", " ", plot_data$visit)
  visit_clean <- gsub("([a-zA-Z])(\\d)", "\\1 \\2", visit_clean)
  plot_data$visit_label <- tools::toTitleCase(visit_clean)

  # Step 4: Preserve visit ordering from pool object
  visit_levels <- unique(plot_data$visit_label)
  plot_data$visit_label <- factor(plot_data$visit_label, levels = visit_levels)

  # Step 5: Build and combine panels
  ...
}
```

### Pattern 3: LSM By-Arm Display with position_dodge
**What:** When displaying LSM estimates (not treatment differences), show two arms per visit using position_dodge.
**When to use:** When user sets `display = "lsm"`.
**Example (verified):**
```r
p_mid <- ggplot(plot_data, aes(x = est, y = visit, xmin = lci, xmax = uci, colour = arm)) +
  geom_linerange(position = position_dodge(width = 0.4), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.4), size = 3, shape = 16) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_colour_manual(values = okabe_ito_colors) +
  scale_y_discrete(limits = rev)
```

### Pattern 4: Dependency Guard (Consistent with efficacy_table)
**What:** Check ggplot2 and patchwork availability before proceeding.
**Example:**
```r
plot_forest <- function(pool_obj, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(
      c("Package {.pkg ggplot2} is required for forest plots.",
        "i" = "Install with {.code install.packages(\"ggplot2\")}."),
      class = c("rbmiUtils_error_dependency", "rbmiUtils_error")
    )
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cli::cli_abort(
      c("Package {.pkg patchwork} is required for forest plots.",
        "i" = "Install with {.code install.packages(\"patchwork\")}."),
      class = c("rbmiUtils_error_dependency", "rbmiUtils_error")
    )
  }
  ...
}
```

### Anti-Patterns to Avoid
- **Using coord_flip() for horizontal orientation:** In ggplot2 4.0+, use `xmin`/`xmax` aesthetics directly with `geom_linerange()`. coord_flip is legacy and causes complications with annotations.
- **Using geom_errorbar/geom_errorbarh for CIs:** These add caps to the whiskers. The decision says "whisker lines without caps" -- use `geom_linerange()` instead.
- **Saving files internally:** The function must return a ggplot2 object, not save a file. Users call `ggsave()` themselves.
- **Moving ggplot2 to Imports:** Keep it in Suggests -- users who only need tables don't need ggplot2 installed.
- **Using geom_pointrange for the CI + point:** While tempting (single geom), it makes it harder to independently control point fill vs. line properties, and complicates the significance highlighting.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Colorblind palette | Custom hex values | `grDevices::palette.colors(palette = "Okabe-Ito")` | Base R built-in since R 4.0; universally recognized standard |
| Multi-panel alignment | grid viewport management | `patchwork::plot_layout(widths = ...)` | Handles alignment, sizing, and ggplot inheritance automatically |
| Visit label cleaning | New regex logic | Copy pattern from `efficacy_table()` | `gsub("_", " ", ...) + gsub("([a-zA-Z])(\\d)", "\\1 \\2", ...) + toTitleCase()` already works |
| P-value formatting | New formatting code | `format_pvalue()` | Already exists in R/formatting.R with tested edge cases |
| Pool object parsing | Manual parameter extraction | `tidy_pool_obj()` | Already handles all parameter naming conventions (ANCOVA, gcomp, underscore visits) |
| Estimate + CI text | sprintf one-off | `format_estimate()` | Already exists in R/formatting.R |

**Key insight:** The bulk of data preparation for `plot_forest()` is already solved by Phase 1's `tidy_pool_obj()` plus the result_helpers (`extract_trt_effects()`, `extract_lsm()`). The new code should focus exclusively on visualization logic.

## Common Pitfalls

### Pitfall 1: Y-axis alignment across patchwork panels
**What goes wrong:** The visit labels on the y-axis don't align between the table panel and the forest plot panel when text widths vary.
**Why it happens:** Each ggplot panel sizes its axes independently.
**How to avoid:** Use `scale_y_discrete(limits = rev)` with identical factor levels in ALL three panels. The `rev` ensures visits appear top-to-bottom. Ensure all panels share the same data source (or at minimum the same factor levels).
**Warning signs:** Visits appearing at different vertical positions in table vs. plot panels.

### Pitfall 2: Patchwork object vs. pure ggplot for user customization
**What goes wrong:** Users expect `plot_forest(obj) + theme(...)` to work on the forest plot specifically.
**Why it happens:** `+ theme()` on a patchwork object applies to the top-level annotation, not sub-plots. Users need `& theme()` to apply to all sub-plots.
**How to avoid:** Document clearly: "Use `& theme()` to customize all panels, or extract individual panels." Alternatively, consider providing `theme` and `text_size` arguments to control common customizations directly.
**Warning signs:** User complaints that theme changes don't take effect.

### Pitfall 3: Visit ordering becomes alphabetical
**What goes wrong:** Visits display as "Week 12, Week 24, Week 4" (alphabetical) instead of pool object order.
**Why it happens:** ggplot2 sorts character vectors alphabetically.
**How to avoid:** Convert visit labels to factors with levels in first-appearance order (same pattern as `efficacy_table()`): `factor(visit_label, levels = unique(visit_label))`.
**Warning signs:** Visits in unexpected order in the plot.

### Pitfall 4: Reference line value mismatch with display mode
**What goes wrong:** Reference line at 0 makes sense for treatment differences but not for LSM display.
**Why it happens:** Zero is the null for treatment differences but meaningless for absolute LSM values.
**How to avoid:** Default reference line to 0 for `display = "trt"`, and `NULL` (no line) for `display = "lsm"`. Allow user override via `ref_value` argument.
**Warning signs:** Confusing reference line placement in LSM plots.

### Pitfall 5: Significance highlighting with non-standard CI levels
**What goes wrong:** Significance is determined by whether CI excludes zero, which depends on the confidence level. If the pool object uses 90% CI, significance markers may disagree with 5% alpha.
**Why it happens:** CI-based significance detection uses whatever confidence level was set in the pool object.
**How to avoid:** Extract `conf.level` from pool object and use it consistently. The CI already reflects the correct level. Document that significance highlighting is based on the CI excluding the reference value, not on p-value thresholds.
**Warning signs:** Inconsistency between p-value column and significance markers.

### Pitfall 6: Testing ggplot2 output in unit tests
**What goes wrong:** Tests try to render plots and check pixel output, leading to fragile tests.
**Why it happens:** Visual output varies across platforms, R versions, and system fonts.
**How to avoid:** Test the ggplot2 object structure (layers, data, aesthetics, scales) rather than rendered output. Use `ggplot2::layer_data()` to extract computed data. Check `class()`, check the data frame underlying the plot, and check that specific geoms are present.
**Warning signs:** Tests failing on CI but passing locally due to font differences.

## Code Examples

### Complete plot_forest() skeleton
```r
# Source: Pattern derived from efficacy_table() + verified ggplot2/patchwork approach
plot_forest <- function(
    pool_obj,
    display = c("trt", "lsm"),
    ref_value = NULL,
    ci_level = NULL,
    arm_labels = NULL,
    title = NULL,
    text_size = 3,
    point_size = 3,
    ...
) {
  # --- Dependency checks ---
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(
      c("Package {.pkg ggplot2} is required for forest plots.",
        "i" = "Install with {.code install.packages(\"ggplot2\")}."),
      class = c("rbmiUtils_error_dependency", "rbmiUtils_error")
    )
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cli::cli_abort(
      c("Package {.pkg patchwork} is required for forest plots.",
        "i" = "Install with {.code install.packages(\"patchwork\")}."),
      class = c("rbmiUtils_error_dependency", "rbmiUtils_error")
    )
  }

  # --- Input validation ---
  if (!inherits(pool_obj, "pool")) {
    cli::cli_abort(
      "Input {.arg pool_obj} must be of class {.cls pool}, not {.cls {class(pool_obj)}}.",
      class = c("rbmiUtils_error_validation", "rbmiUtils_error")
    )
  }

  display <- match.arg(display)

  # --- Metadata ---
  if (is.null(ci_level)) {
    ci_level <- pool_obj$conf.level %||% 0.95
  }
  if (is.null(ref_value)) {
    ref_value <- if (display == "trt") 0 else NULL
  }

  # --- Data preparation ---
  tidy_df <- tidy_pool_obj(pool_obj)

  # Filter by display mode
  if (display == "trt") {
    plot_data <- tidy_df[tidy_df$parameter_type == "trt", ]
  } else {
    plot_data <- tidy_df[tidy_df$parameter_type == "lsm", ]
  }

  # Handle NA visits
  plot_data <- plot_data[!is.na(plot_data$visit), ]

  # Clean visit labels (reuse efficacy_table pattern)
  visit_clean <- gsub("_", " ", plot_data$visit)
  visit_clean <- gsub("([a-zA-Z])(\\d)", "\\1 \\2", visit_clean)
  plot_data$visit_label <- factor(
    tools::toTitleCase(visit_clean),
    levels = unique(tools::toTitleCase(visit_clean))
  )

  # Significance flag (CI excludes reference value)
  if (!is.null(ref_value)) {
    plot_data$significant <- (plot_data$lci > ref_value) | (plot_data$uci < ref_value)
  } else {
    plot_data$significant <- FALSE
  }

  # Format text columns
  plot_data$est_ci_label <- format_estimate(plot_data$est, plot_data$lci, plot_data$uci)
  plot_data$pval_label <- ifelse(
    is.na(plot_data$pval), "\u2014",
    format_pvalue(plot_data$pval)
  )

  # --- Build panels ---
  # ... (panel construction code)

  # --- Combine with patchwork ---
  combined <- p_left + p_mid + p_right +
    patchwork::plot_layout(widths = c(3, 4, 1.5))

  if (!is.null(title)) {
    combined <- combined + patchwork::plot_annotation(title = title)
  }

  combined
}
```

### Significance highlighting approach (recommended)
```r
# Use alpha (transparency) to highlight: significant = full opacity, non-significant = faded
# This is subtle but effective and doesn't change shape (per decision: uniform filled circles)
ggplot2::geom_point(
  ggplot2::aes(alpha = significant),
  size = point_size, shape = 16
) +
ggplot2::scale_alpha_manual(
  values = c("TRUE" = 1.0, "FALSE" = 0.4),
  guide = "none"
)
```

Alternative: use filled vs open circles (shape 16 vs 21):
```r
ggplot2::geom_point(
  ggplot2::aes(shape = significant),
  size = point_size
) +
ggplot2::scale_shape_manual(
  values = c("TRUE" = 16, "FALSE" = 1),  # 16 = filled, 1 = open
  guide = "none"
)
```

### Okabe-Ito palette extraction (no external dependency)
```r
# Built into base R since 4.0 -- no package needed
okabe_ito <- grDevices::palette.colors(palette = "Okabe-Ito")
# Named vector: black, orange, skyblue, bluegreen, yellow, blue, vermillion, reddishpurple, grey

# For two-arm clinical trial: use blue (#0072B2) and vermillion (#D55E00)
# These are maximally distinguishable for all forms of color vision deficiency
forest_colors <- c(
  "ref" = unname(okabe_ito[6]),    # blue: #0072B2
  "alt" = unname(okabe_ito[7])     # vermillion: #D55E00
)
```

### Clinical theme (white background, minimal gridlines)
```r
theme_forest <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}
```

### Testing ggplot2 objects (verified pattern)
```r
test_that("plot_forest returns patchwork/ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  mock_pool <- make_mock_pool()
  p <- plot_forest(mock_pool)

  # Check class hierarchy
  expect_s3_class(p, "patchwork")
  expect_s3_class(p, "ggplot")  # patchwork inherits from ggplot
})

test_that("plot_forest contains expected layers", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  mock_pool <- make_mock_pool()
  p <- plot_forest(mock_pool)

  # Extract the forest plot panel (middle panel)
  # Patchwork stores patches in p$patches$plots list
  # The first plot added (left panel) is in patches, main plot is p itself
  # Test data structure instead of visual output
})

test_that("plot data contains correct visits", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  mock_pool <- make_mock_pool()
  p <- plot_forest(mock_pool)

  # Access data from the patchwork object
  # For structural testing, better to test the internal data prep function separately
})
```

### Mock pool object (reuse from efficacy_table tests)
```r
# Same make_mock_pool() pattern from test-efficacy_table.R
make_mock_pool <- function() {
  pool_obj <- list(
    pars = list(
      trt_Week4 = list(est = -2.5, se = 0.8, ci = c(-4.1, -0.9), pvalue = 0.002),
      lsm_ref_Week4 = list(est = 10.0, se = 0.5, ci = c(9.0, 11.0), pvalue = NA),
      lsm_alt_Week4 = list(est = 7.5, se = 0.6, ci = c(6.3, 8.7), pvalue = NA),
      trt_Week8 = list(est = -1.0, se = 1.2, ci = c(-3.4, 1.4), pvalue = 0.42)
    ),
    conf.level = 0.95,
    alternative = "two.sided",
    N = 100,
    method = "rubin"
  )
  class(pool_obj) <- "pool"
  pool_obj
}
```

## Recommended Decisions (Claude's Discretion Items)

### Plot orientation: Horizontal
**Decision:** Horizontal forest plot (visits on y-axis, estimate on x-axis).
**Rationale:** Standard clinical trial convention (CSR, FDA submissions). Vertical forest plots are unconventional and would require users to rotate their heads to read visit labels. The KHstats, tern, forestploter, and ggforestplot packages all default to horizontal.

### Display mode default: Treatment difference
**Decision:** Default `display = "trt"` showing treatment differences only. User can switch to `display = "lsm"` for LSM by arm.
**Rationale:** Per CONTEXT.md decisions: "Default to treatment difference when user doesn't specify." Treatment differences are the primary regulatory interest.

### Visit ordering: First-appearance from pool object
**Decision:** Use first-appearance order from `tidy_pool_obj()` output, converted to factor levels.
**Rationale:** Matches `efficacy_table()` behavior (already tested and verified to preserve pool object order, not alphabetical). Same pattern: `factor(visit_label, levels = unique(visit_label))`.

### Multi-endpoint: Single endpoint per call (no faceting)
**Decision:** `plot_forest()` handles a single pool object (single endpoint). No built-in faceting.
**Rationale:** Pool objects from rbmi are per-endpoint. Multi-endpoint support would require `combine_results()` output, which is a different input type. Users can use patchwork themselves: `plot_forest(pool1) / plot_forest(pool2)`.

### Table panel columns: Estimate (CI) and P-value
**Decision:** Left panel shows visit label + "Estimate (LCI, UCI)" text. Right panel shows p-value.
**Rationale:** This is the standard clinical forest plot layout. Matches what CSR Table 14.2.x pairs with.

### Table panel implementation: Patchwork three-panel
**Decision:** Use patchwork with `plot_layout(widths = c(3, 4, 1.5))` combining three ggplot objects.
**Rationale:** Verified working on installed versions. Clean alignment. Returns object that inherits from ggplot (user can customize). The pure-ggplot `clip = "off"` approach is fragile with varying text widths.

### Significance highlighting: Filled vs. open circles
**Decision:** Significant results (CI excludes reference value) use filled circles (shape 16). Non-significant use open circles (shape 1).
**Rationale:** Clear visual distinction without changing color. Works in both color and B&W printing. Subtle but effective. The alternative (alpha/transparency) can look like a rendering artifact.

### API arguments
**Decision:** Moderate surface -- essential arguments only:
```r
plot_forest(
  pool_obj,              # Required: rbmi pool object
  display = "trt",       # "trt" or "lsm"
  ref_value = NULL,      # Reference line value (default: 0 for trt, NULL for lsm)
  ci_level = NULL,       # CI level label (extracted from pool_obj if NULL)
  arm_labels = NULL,     # Named c(ref = "...", alt = "...") for LSM mode
  title = NULL,          # Plot title
  text_size = 3,         # geom_text size for table panels
  point_size = 3         # geom_point size
)
```
**Rationale:** Mirrors `efficacy_table()` argument style. Users who need more control can use `& theme()` on the returned patchwork object.

### Color palette: Okabe-Ito (blue + vermillion)
**Decision:** Use Okabe-Ito blue (#0072B2) for reference arm and vermillion (#D55E00) for treatment arm.
**Rationale:** Okabe-Ito is the gold standard for colorblind-friendly scientific palettes. Blue and vermillion (orange-red) are the most distinguishable pair across all forms of color vision deficiency. Available via `grDevices::palette.colors()` in base R -- no external dependency.

### Point size and line weight defaults
**Decision:** Point size = 3, line width = 0.6 for CI whiskers, 0.5 for reference line.
**Rationale:** Standard ggplot2 defaults are slightly too small for publication. Size 3 is legible at standard CSR dimensions (typically 6-8 inches wide).

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `coord_flip()` for horizontal plots | Direct `xmin`/`xmax` mapping | ggplot2 3.3+ (2020) | Simpler code, better annotation support |
| `geom_errorbarh()` | `geom_errorbar(orientation = "y")` | ggplot2 3.3+ | `geom_errorbarh` deprecated |
| gridExtra/cowplot for multi-panel | patchwork | patchwork 1.0 (2019), now 1.3.2 | Better API, `& theme()` operator, gt table support |
| Custom colorblind hex codes | `grDevices::palette.colors("Okabe-Ito")` | R 4.0 (2020) | Built-in, no dependency |
| S3 ggplot class checks | `ggplot2::is_ggplot()` | ggplot2 3.5.2 | `is.ggplot()` deprecated |

**Deprecated/outdated:**
- `geom_errorbarh()`: Use `geom_linerange()` or `geom_errorbar(orientation = "y")` instead
- `is.ggplot()`: Use `ggplot2::is_ggplot()` (changed in ggplot2 3.5.2)
- `cowplot::plot_grid()` for simple compositions: patchwork's `+` operator is simpler and better supported

## Open Questions

1. **Patchwork internal structure for testing**
   - What we know: Patchwork objects have `$patches$plots` list containing sub-plots. The "main" plot is the patchwork object itself.
   - What's unclear: The exact internal structure for extracting individual panel data in tests.
   - Recommendation: Test data preparation logic separately (extract the data prep into a testable internal function `prepare_forest_data()`), and test the composed plot at the structural level (class checks, layer counts).

2. **Column header alignment in table panels**
   - What we know: Table panels use `geom_text()` + `theme_void()`. Column headers (e.g., "Estimate (95% CI)") can be set via `labs(title = ...)` on each sub-plot.
   - What's unclear: Whether patchwork perfectly aligns these titles with the forest plot header.
   - Recommendation: Test visually during implementation. If alignment is poor, use `patchwork::plot_annotation()` for the overall title and skip per-panel titles.

3. **LSM arm label extraction from pool object**
   - What we know: `tidy_pool_obj()` extracts `lsm_type` as "ref" or "alt" for standard ANCOVA, or actual arm names for gcomp.
   - What's unclear: Whether `arm_labels` argument should map to `lsm_type` values or to actual arm names from the pool object.
   - Recommendation: Default labels from `lsm_type` ("Reference"/"Treatment"), with `arm_labels = c(ref = "Placebo", alt = "Drug A")` override. Same pattern as `efficacy_table()`.

## Sources

### Primary (HIGH confidence)
- ggplot2 4.0.1 installed locally -- verified `geom_linerange()`, `geom_point()`, horizontal orientation
- patchwork 1.3.2 installed locally -- verified `plot_layout(widths = ...)`, `& theme()` operator
- `grDevices::palette.colors("Okabe-Ito")` verified locally -- 9 colors available in base R
- Existing codebase: `tidy_pool_obj()`, `efficacy_table()`, `format_pvalue()`, `format_estimate()` patterns verified
- [ggplot2 interval geom reference](https://ggplot2.tidyverse.org/reference/geom_linerange.html) -- `xmin`/`xmax` semantics confirmed

### Secondary (MEDIUM confidence)
- [KHstats Annotated Forest Plots](https://www.khstats.com/blog/forest-plots/) -- three-panel patchwork approach confirmed
- [patchwork layout documentation](https://patchwork.data-imaginist.com/articles/guides/layout.html) -- width control syntax confirmed
- [Okabe-Ito palette reference](https://easystats.github.io/see/reference/scale_color_okabeito.html) -- hex values and colorblind properties confirmed
- [tern g_forest](https://insightsengineering.github.io/tern/latest-tag/reference/g_forest.html) -- industrial pharma forest plot pattern reference

### Tertiary (LOW confidence)
- None -- all findings verified with primary or secondary sources

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- ggplot2 and patchwork verified locally, both installed
- Architecture: HIGH -- three-panel patchwork pattern verified with working code
- Pitfalls: HIGH -- derived from practical testing and codebase patterns
- API design: HIGH -- consistent with existing efficacy_table() pattern in codebase
- Color palette: HIGH -- Okabe-Ito verified from base R

**Research date:** 2026-02-08
**Valid until:** 2026-05-08 (ggplot2 4.x stable; patchwork 1.x stable)
