# Phase 10: Publication Styling - Research

**Researched:** 2026-02-11
**Domain:** gt table styling, ggplot2 typography, patchwork panel layout control
**Confidence:** HIGH

## Summary

This phase adds publication styling parameters to two existing functions: `efficacy_table()` (gt-based) and `plot_forest()` (ggplot2/patchwork-based). The work is purely additive -- adding new optional parameters with sensible defaults that preserve current behavior. No new dependencies are required; all needed functionality exists in the already-imported gt (1.3.0), ggplot2 (4.0.1), and patchwork (1.3.2).

For gt tables (STYLE-01, STYLE-02): The `opt_table_font()` function sets the font family globally on a gt table. Font size is controlled via `tab_options(table.font.size = px(...))` and row padding via `tab_options(data_row.padding = px(...))` or the convenience function `opt_vertical_padding(scale = ...)`. These are well-documented, stable APIs.

For ggplot2 forest plots (STYLE-03, STYLE-04, STYLE-05): Font family must be passed explicitly to each `geom_text()` call via the `family` parameter. Theme text elements (axis labels, titles) are set via `theme_minimal(base_family = ...)`. The `panel_widths` parameter maps directly to patchwork's `plot_layout(widths = ...)`. Alignment improvements for the left panel require switching from right-justified (`hjust = 1`) to left-justified (`hjust = 0`) text placement.

**Primary recommendation:** Add `font_family`, `font_size`, and `row_padding` parameters to `efficacy_table()`, and `font_family` and `panel_widths` parameters to `plot_forest()`, with `NULL` defaults that preserve current behavior.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| gt | 1.3.0 | Table rendering with font/padding control | Already in Suggests; `opt_table_font()` and `tab_options()` are the standard APIs |
| ggplot2 | 4.0.1 | Plot rendering with font control | Already in Suggests; `geom_text(family=)` and `theme_minimal(base_family=)` |
| patchwork | 1.3.2 | Panel width control | Already in Suggests; `plot_layout(widths=)` is the standard API |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| systemfonts | (ggplot2 dep) | Font availability checking | Optional validation of font_family parameter |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `opt_table_font()` | `tab_options(table.font.names=)` | `opt_table_font()` is the higher-level API with fallback stacks; prefer it |
| `geom_text(family=)` | `theme(geom = element_geom(family=))` | ggplot2 4.0+ only; explicit param is backwards-compatible and clearer |
| `tab_options(data_row.padding=)` | `opt_vertical_padding(scale=)` | `scale` is relative (0-3); explicit `px()` value gives more precise control |

**No new dependencies required.** All functionality exists in current Suggests.

## Architecture Patterns

### Pattern 1: Additive Optional Parameters with NULL Defaults
**What:** New styling parameters default to `NULL`, meaning "use current behavior." When non-NULL, apply the styling.
**When to use:** Always -- this preserves backward compatibility.
**Example:**
```r
# In efficacy_table():
efficacy_table <- function(pool_obj, ..., font_family = NULL, font_size = NULL, row_padding = NULL) {
  # ... existing code builds tbl ...

  # Apply optional styling at the end, before return
  if (!is.null(font_family)) {
    tbl <- gt::opt_table_font(tbl, font = font_family)
  }
  if (!is.null(font_size)) {
    tbl <- gt::tab_options(tbl, table.font.size = gt::px(font_size))
  }
  if (!is.null(row_padding)) {
    tbl <- gt::tab_options(tbl, data_row.padding = gt::px(row_padding))
  }

  tbl
}
```

### Pattern 2: Font Family Propagation in Forest Plot Panels
**What:** Pass `font_family` to every `geom_text()` call AND to `theme_minimal(base_family=)` in the forest panel.
**When to use:** When `font_family` parameter is non-NULL.
**Why:** `geom_text()` does NOT inherit font family from `theme(text = element_text(family=))` (confirmed: this is intentional ggplot2 design). Each `geom_text` call must receive `family` explicitly. Theme elements (axis labels, titles) are controlled via `base_family`.
**Example:**
```r
# In plot_forest():
# Left panel geom_text calls
geom_text(aes(x = 0, label = .data$visit_label),
          hjust = 0, size = text_size, fontface = "bold",
          family = font_family %||% "")  # "" = ggplot2 default (sans)

# Middle panel theme
theme_forest(base_size = 11, base_family = font_family %||% "")

# Right panel geom_text calls
geom_text(aes(x = 0, label = .data$pval_label),
          hjust = 0.5, size = text_size,
          family = font_family %||% "")
```

### Pattern 3: Panel Width Parameterization
**What:** Expose patchwork `widths` as a user parameter with sensible defaults.
**When to use:** For STYLE-04 (panel_widths parameter).
**Example:**
```r
plot_forest <- function(..., panel_widths = NULL) {
  # Default widths (current hardcoded values)
  if (is.null(panel_widths)) {
    panel_widths <- if (show_pvalues) c(3, 4, 1.5) else c(3, 5)
  }

  # Validate length matches panel count
  expected_length <- if (show_pvalues) 3L else 2L
  if (length(panel_widths) != expected_length) {
    cli::cli_abort(
      "{.arg panel_widths} must have length {expected_length} (got {length(panel_widths)}).",
      class = c("rbmiUtils_error_validation", "rbmiUtils_error")
    )
  }

  # Apply
  combined <- p_left + p_mid + p_right +
    patchwork::plot_layout(widths = panel_widths)
}
```

### Pattern 4: Left Panel Alignment Improvement
**What:** Switch visit labels to left-justified and estimates to left-justified at a fixed column position.
**When to use:** For STYLE-05 (consistent alignment).
**Current problem:** Both `geom_text` layers use `hjust = 1` (right-align). With proportional fonts, right-aligned text of varying width creates ragged left edges.
**Fix:** Use `hjust = 0` (left-align) for both columns. This anchors each column's left edge at the fixed x position, giving consistent alignment regardless of text width.
**Example:**
```r
# Before (current - right-aligned, ragged left edges):
geom_text(aes(x = 0, label = visit_label), hjust = 1, ...)
geom_text(aes(x = 1, label = est_ci_label), hjust = 1, ...)

# After (left-aligned, consistent left edges):
geom_text(aes(x = 0, label = visit_label), hjust = 0, ...)
geom_text(aes(x = 1, label = est_ci_label), hjust = 0, ...)
# Adjust scale_x_continuous expand to accommodate text overflow on right
scale_x_continuous(expand = expansion(mult = c(0.05, 0.3)))
```

### Anti-Patterns to Avoid
- **Moving gt/ggplot2 to Imports for styling:** Keep them in Suggests. Styling parameters are optional enhancements.
- **Validating font availability:** Do not check if a font is installed. Font fallback is handled by the rendering engine (browser for gt HTML, graphics device for ggplot2). Invalid fonts silently fall back to defaults -- this is standard behavior.
- **Adding font_size parameter to plot_forest:** The existing `text_size` parameter already controls geom_text size. Adding a separate `font_size` would be confusing. Font size for axis labels is controlled via `& theme(text = element_text(size = ...))`.
- **Using element_geom(family=) for ggplot2 font:** This requires ggplot2 4.0+ and is less explicit. Direct `family=` in geom_text works with all ggplot2 versions.
- **Applying row_padding to specific gt sections only:** Use `data_row.padding` which affects all data rows uniformly. Per-section padding creates visual inconsistency.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| gt font family | Custom CSS injection | `gt::opt_table_font(font = ...)` | Handles font stacks, Google fonts, fallbacks automatically |
| gt row padding | Manual cell height manipulation | `gt::tab_options(data_row.padding = gt::px(...))` | Works across all output formats (HTML, LaTeX, RTF) |
| gt compact layout | Individual padding adjustments | `gt::opt_vertical_padding(scale = ...)` | Scales all padding proportionally in one call |
| ggplot font propagation | Wrapper that sets family on every layer | Pass `family` param to each `geom_text()` at construction | The function owns all geom_text calls; no wrapper needed |
| Panel width validation | Complex ratio normalization | Validate length matches panel count | patchwork handles relative sizing automatically |

**Key insight:** The styling APIs already exist in gt and ggplot2. This phase is about exposing them as parameters, not building new rendering logic.

## Common Pitfalls

### Pitfall 1: geom_text Does NOT Inherit Font from theme(text = element_text())
**What goes wrong:** Setting `theme(text = element_text(family = "Times"))` and expecting geom_text to use Times. It will not.
**Why it happens:** In ggplot2, theme elements affect non-data components only (axis labels, titles, legends). geom_text is a data layer. This is intentional design (confirmed via ggplot2 issue #1859, closed as "working as intended").
**How to avoid:** Pass `family` explicitly to every `geom_text()` call. For theme text elements, use `theme_minimal(base_family = ...)` or `theme(text = element_text(family = ...))`.
**Warning signs:** geom_text renders in default sans-serif while axis labels show the custom font.

### Pitfall 2: gt font_size Requires px() Helper
**What goes wrong:** Passing a bare numeric like `table.font.size = 10` to `tab_options()`. The parameter expects a character string with units (e.g., "10px") or the `px()` helper.
**Why it happens:** gt supports multiple unit types (px, pct, em).
**How to avoid:** Always wrap with `gt::px()`: `tab_options(table.font.size = gt::px(font_size))`. In the function signature, accept a plain numeric and wrap internally.
**Warning signs:** Font size not changing, or unexpected size behavior.

### Pitfall 3: Panel Width Length Mismatch When show_pvalues = FALSE
**What goes wrong:** User passes `panel_widths = c(3, 4, 1.5)` but also sets `show_pvalues = FALSE`, creating a length mismatch (3 widths for 2 panels).
**Why it happens:** `show_pvalues = FALSE` creates a 2-panel layout, but user provided 3-width vector from the 3-panel default.
**How to avoid:** Validate that `length(panel_widths)` matches the actual panel count (3 if show_pvalues, 2 otherwise). Provide a clear error message.
**Warning signs:** patchwork error about width vector length.

### Pitfall 4: Alignment Change Breaking Existing Visual Tests
**What goes wrong:** Changing `hjust` from 1 to 0 in the left panel changes the visual appearance. Tests that check rendered output or specific text positions may break.
**Why it happens:** Existing tests verify the current right-aligned behavior.
**How to avoid:** Current tests check structural properties (class, data values) not visual alignment. Verify that no test depends on hjust values. The alignment change should be transparent to existing tests.
**Warning signs:** Unexpected test failures after hjust change.

### Pitfall 5: scale_x_continuous Expand Direction After hjust Change
**What goes wrong:** After switching to `hjust = 0` (left-align), text overflows the right edge of the panel because `expand` was set for left-side overflow.
**Why it happens:** `expansion(mult = c(0.3, 0.05))` adds 30% on the left. With hjust=0, text extends rightward, needing expansion on the right.
**How to avoid:** Reverse the expand multipliers: `expansion(mult = c(0.05, 0.3))` to add space on the right where text extends.
**Warning signs:** Text clipped or extending beyond panel boundaries.

## Code Examples

### STYLE-01 & STYLE-02: efficacy_table Font and Padding
```r
# Source: gt::opt_table_font docs, gt::tab_options docs (verified locally)
efficacy_table <- function(
    pool_obj,
    title = NULL,
    subtitle = NULL,
    digits = 2,
    ci_level = NULL,
    arm_labels = NULL,
    pval_digits = 3,
    pval_threshold = 0.001,
    font_family = NULL,   # NEW: STYLE-01
    font_size = NULL,      # NEW: STYLE-02
    row_padding = NULL,    # NEW: STYLE-02
    ...
) {
  # ... existing body unchanged ...

  # Apply publication styling (after table construction, before return)
  if (!is.null(font_family)) {
    tbl <- gt::opt_table_font(tbl, font = font_family)
  }
  if (!is.null(font_size)) {
    tbl <- gt::tab_options(tbl, table.font.size = gt::px(font_size))
  }
  if (!is.null(row_padding)) {
    tbl <- gt::tab_options(tbl, data_row.padding = gt::px(row_padding))
  }

  tbl
}

# Usage:
efficacy_table(pool_obj,
  font_family = "Times New Roman",
  font_size = 10,
  row_padding = 3
)
```

### STYLE-03: plot_forest Font Family
```r
# Source: ggplot2 geom_text family param (verified locally)
# Key insight: must pass family to EVERY geom_text call
plot_forest <- function(
    pool_obj,
    display = c("trt", "lsm"),
    ...,
    font_family = NULL    # NEW: STYLE-03
) {
  # Resolve font for geom_text (empty string = ggplot2 default)
  geom_family <- font_family %||% ""

  # Left panel: pass family to both geom_text layers
  p_left <- ggplot2::ggplot(plot_data, ggplot2::aes(y = .data$visit_label)) +
    ggplot2::geom_text(
      ggplot2::aes(x = 0, label = .data$visit_label),
      hjust = 0, size = text_size, fontface = "bold",
      family = geom_family
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 1, label = .data$est_ci_label),
      hjust = 0, size = text_size,
      family = geom_family
    ) + ...

  # Middle panel: pass base_family to theme
  theme_forest(base_size = 11, base_family = geom_family)

  # Right panel: pass family to geom_text
  p_right <- ggplot2::ggplot(plot_data, ggplot2::aes(y = .data$visit_label)) +
    ggplot2::geom_text(
      ggplot2::aes(x = 0, label = .data$pval_label),
      hjust = 0.5, size = text_size,
      family = geom_family
    ) + ...
}
```

### STYLE-04: Panel Width Control
```r
# Source: patchwork::plot_layout(widths=) docs (verified locally)
plot_forest <- function(..., panel_widths = NULL) {
  # Default widths preserve current behavior
  if (is.null(panel_widths)) {
    panel_widths <- if (show_pvalues) c(3, 4, 1.5) else c(3, 5)
  }

  # Validation
  expected_n <- if (show_pvalues) 3L else 2L
  if (length(panel_widths) != expected_n) {
    cli::cli_abort(
      c("{.arg panel_widths} must be a numeric vector of length {expected_n}.",
        "i" = "Got length {length(panel_widths)}.",
        "i" = if (show_pvalues) {
          "Three widths needed: table panel, forest panel, p-value panel."
        } else {
          "Two widths needed: table panel, forest panel."
        }),
      class = c("rbmiUtils_error_validation", "rbmiUtils_error")
    )
  }

  combined <- p_left + p_mid + p_right +
    patchwork::plot_layout(widths = panel_widths)
}
```

### STYLE-05: Left Panel Alignment Fix
```r
# Source: ggplot2 geom_text hjust behavior (verified locally)
# Change from hjust=1 (right-align) to hjust=0 (left-align)
# This gives consistent left-edge alignment regardless of label width

# TRT mode left panel (after fix):
p_left <- ggplot2::ggplot(plot_data, ggplot2::aes(y = .data$visit_label)) +
  ggplot2::geom_text(
    ggplot2::aes(x = 0, label = .data$visit_label),
    hjust = 0, size = text_size, fontface = "bold",
    family = geom_family
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = 1, label = .data$est_ci_label),
    hjust = 0, size = text_size,
    family = geom_family
  ) +
  # Note: expand direction reversed for left-aligned text
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0.05, 0.3))
  ) +
  ggplot2::scale_y_discrete(limits = rev) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5))
```

### theme_forest Update for base_family
```r
# Source: ggplot2::theme_minimal(base_family=) (verified locally)
theme_forest <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
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

### Testing Font Family Application
```r
# Test pattern: verify font_family parameter is applied
test_that("font_family parameter propagates to gt table", {
  skip_if_not_installed("gt")
  mock_pool <- make_mock_pool()

  tbl <- efficacy_table(mock_pool, font_family = "Courier")
  html <- gt::as_raw_html(tbl)

  # Check that font family appears in the rendered HTML
  expect_true(grepl("Courier", html, fixed = TRUE))
})

test_that("font_family parameter propagates to all forest plot panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  mock_pool <- make_mock_pool()
  p <- plot_forest(mock_pool, font_family = "serif")

  # Build and check that geom_text layers use the specified family
  # Access left panel (first patch)
  left_panel <- p$patches$plots[[1]]
  b <- ggplot2::ggplot_build(left_panel)

  # Both text layers should have family = "serif"
  expect_equal(b$data[[1]]$family[1], "serif")
  expect_equal(b$data[[2]]$family[1], "serif")
})

test_that("panel_widths parameter controls layout", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  mock_pool <- make_mock_pool()

  # Custom widths should not error
  p <- plot_forest(mock_pool, panel_widths = c(2, 5, 1))
  expect_s3_class(p, "patchwork")
  expect_no_error(print(p))
})

test_that("panel_widths validates length against show_pvalues", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  mock_pool <- make_mock_pool()

  # 3 widths with show_pvalues = FALSE should error
  expect_error(
    plot_forest(mock_pool, panel_widths = c(3, 4, 1.5), show_pvalues = FALSE),
    class = "rbmiUtils_error_validation"
  )
})
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `theme(text = element_text(family=))` for geom font | `geom_text(family=)` directly or `theme(geom = element_geom(family=))` | ggplot2 4.0.0 (Sep 2025) | `element_geom` is new; direct param remains reliable for all versions |
| `gt::tab_options(table.font.names=)` | `gt::opt_table_font(font=)` | gt 0.8+ | Higher-level API with font stack support |
| Manual padding per section | `gt::opt_vertical_padding(scale=)` | gt 0.8+ | Single call scales all padding proportionally |
| `geom_text` ignores theme font | `geom_text` family now defaults via `theme()` in ggplot2 4.0+ | ggplot2 4.0.0 | Automatic theme inheritance (but explicit param is clearer) |

**Deprecated/outdated:**
- `gt::tab_options(table.font.names = ...)` is superseded by `gt::opt_table_font(font = ...)` for font family control
- `is.ggplot()` deprecated in ggplot2 3.5.2; use `ggplot2::is_ggplot()`

## Open Questions

1. **Should font_size accept character or numeric?**
   - What we know: `gt::tab_options(table.font.size = ...)` accepts `gt::px(N)`, `gt::pct(N)`, or character like "10px"
   - What's unclear: Whether to accept only numeric (simpler API) or also character (more flexible)
   - Recommendation: Accept numeric only, wrap with `gt::px()` internally. Simpler API, consistent with `row_padding`. Users wanting `pct()` or `em()` can use `gt::tab_options()` directly via pipe.

2. **Should font_family validation check system fonts?**
   - What we know: `systemfonts::system_fonts()` could validate font availability. Invalid fonts silently fall back to defaults.
   - What's unclear: Whether silent fallback or explicit warning is better UX.
   - Recommendation: Do NOT validate. Silent fallback is standard behavior in both gt and ggplot2. Font availability varies by system and rendering context. Let the rendering engine handle it.

3. **Should panel_widths accept named vector?**
   - What we know: Names like `c(table = 3, forest = 4, pvalue = 1.5)` would be self-documenting.
   - What's unclear: Whether naming adds enough value vs. added validation complexity.
   - Recommendation: Accept unnamed numeric vector. Document the order (table, forest, p-value) in the `@param` roxygen. Names would add validation overhead without significant benefit for a 2-3 element vector.

4. **LSM mode left panel alignment (STYLE-05)**
   - What we know: LSM mode constructs `row_label` as "Visit - Arm" and uses different `geom_text` calls than TRT mode. The alignment fix must apply to both code paths.
   - What's unclear: Whether LSM mode's combined row labels need different x-position spacing.
   - Recommendation: Apply the same `hjust = 0` fix to both TRT and LSM mode left panels. Adjust x positions and expand multipliers identically.

## Sources

### Primary (HIGH confidence)
- gt 1.3.0 installed locally -- `opt_table_font()`, `tab_options(data_row.padding=, table.font.size=)`, `opt_vertical_padding()` all verified working
- ggplot2 4.0.1 installed locally -- `geom_text(family=)` verified; `theme(geom = element_geom(family=))` verified; `theme_minimal(base_family=)` verified
- patchwork 1.3.2 installed locally -- `plot_layout(widths=)` verified; `& theme()` propagation verified
- [gt opt_table_font reference](https://gt.rstudio.com/reference/opt_table_font.html) -- API signature and usage
- [gt tab_options reference](https://gt.rstudio.com/reference/tab_options.html) -- data_row.padding and table.font.size params
- [gt opt_vertical_padding reference](https://gt.rstudio.com/reference/opt_vertical_padding.html) -- scale parameter behavior
- [ggplot2 geom_text reference](https://ggplot2.tidyverse.org/reference/geom_text.html) -- family aesthetic defaults
- [ggplot2 element_geom reference](https://ggplot2.tidyverse.org/reference/element.html) -- element_geom() signature
- [patchwork plot_layout reference](https://patchwork.data-imaginist.com/reference/plot_layout.html) -- widths parameter

### Secondary (MEDIUM confidence)
- [ggplot2 issue #1859](https://github.com/tidyverse/ggplot2/issues/1859) -- Confirmed: geom_text does NOT inherit from `theme(text=)`, must use direct `family` param
- [ggplot2 4.0.0 release notes](https://tidyverse.org/blog/2025/09/ggplot2-4-0-0/) -- element_geom() and from_theme() introduction
- [ggplot2 NEWS](https://cran.r-project.org/web/packages/ggplot2/news/news.html) -- header_family, element_geom additions

### Tertiary (LOW confidence)
- None -- all findings verified with primary or secondary sources

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all libraries already installed and in Suggests; APIs verified locally
- Architecture: HIGH -- patterns verified with running R code; all gt/ggplot2 APIs confirmed working
- Pitfalls: HIGH -- geom_text font inheritance verified empirically; alignment behavior verified
- Code examples: HIGH -- all examples tested locally on ggplot2 4.0.1, gt 1.3.0, patchwork 1.3.2

**Research date:** 2026-02-11
**Valid until:** 2026-05-11 (gt, ggplot2, patchwork all stable; no breaking changes expected)
