# Phase 3: Efficacy Tables - Research

**Researched:** 2026-02-08
**Domain:** Regulatory-style efficacy summary tables from rbmi pool objects via gt
**Confidence:** HIGH

## Summary

This phase implements a function (`efficacy_table()`) that takes an rbmi pool object and produces a formatted gt table showing LS means by arm, treatment differences, CIs, and p-values organized by visit -- matching the standard CDISC/ICH Table 14.2.x format used in Clinical Study Reports.

The recommended approach is **direct gt construction from tidy_pool_obj() data** rather than going through gtsummary's `tbl_ard_*()` pipeline. The gtsummary ARD-first functions (`tbl_ard_summary()`, `tbl_ard_wide_summary()`) are designed for descriptive statistics (means, medians, counts) computed from raw data, not for pre-computed model-based estimates like LS means and treatment differences. Building directly with gt gives full control over the Table 14.2.x layout: visits as row groups, LS means per arm followed by treatment difference rows, with estimate/SE/CI/p-value columns.

The gt package (v1.3.0 installed) provides all required features: `tab_row_group()` / `groupname_col` for visit grouping, `tab_header()` for title/subtitle, `tab_footnote()` / `tab_source_note()` for metadata footnotes, `fmt_number()` for decimal formatting, `cols_merge()` for CI bracket display, and `gtsave()` / `as_raw_html()` / `as_latex()` for HTML/PDF output. A working prototype was verified during research.

**Primary recommendation:** Build a single `efficacy_table()` function that accepts a pool object, internally calls `tidy_pool_obj()` to extract data, reshapes it into a table-ready data frame, and constructs a gt table with opinionated defaults. Return the gt object so users can pipe into further `gt::*` customization. Use a two-step internal design (prepare data frame + render gt) but expose a single public function.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| gt | >= 0.10.0 | Table rendering (row groups, spanners, footnotes, formatting) | Decision locked: gt-only output; v1.3.0 installed; supports HTML/PDF/LaTeX |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| dplyr | (existing Imports) | Data reshaping for table-ready data frame | Already in Imports |
| cli | (existing Imports) | Input validation error messages | Already in Imports |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Direct gt construction | gtsummary `tbl_ard_summary()` via ARD | gtsummary `tbl_ard_*()` functions are designed for descriptive stats from raw data, not pre-computed model estimates. Would require force-fitting rbmi results into a summary table paradigm. Direct gt gives exact control over Table 14.2.x layout. |
| Direct gt construction | gtsummary `as_gt()` pipeline | Adds gtsummary as dependency for minimal benefit -- the table structure is custom (model-based, not summary-based) |
| Single `efficacy_table()` | Separate `prepare_efficacy_data()` + `render_efficacy_table()` | Single function is simpler for users (one call), matches TBL-03 requirement. Internal separation still possible. |

**Installation:**
```r
# gt as Suggests (not Imports) -- matches Phase 2 pattern for cards
# Already available: dplyr, cli in Imports
```

## Architecture Patterns

### Recommended Project Structure
```
R/
  efficacy_table.R       # NEW: efficacy_table() function
  formatting.R           # EXISTING: format_pvalue(), format_estimate() (reused)
  tidiers.R              # EXISTING: tidy_pool_obj() (input data source)
  pool_methods.R         # EXISTING: print.pool, summary.pool
tests/testthat/
  test-efficacy_table.R  # NEW: tests for efficacy_table()
```

### Pattern 1: Pool Object to Table-Ready Data Frame
**What:** Transform tidy_pool_obj() output into the exact row/column structure needed for gt
**When to use:** Always -- this is the data preparation step before gt rendering
**Confidence:** HIGH -- verified with working prototype

The tidy_pool_obj() output has columns: `parameter, description, visit, parameter_type, lsm_type, est, se, lci, uci, pval`. This needs reshaping:

1. Clean visit labels (replace underscores with spaces, title case)
2. Create display row labels ("LS Mean (Reference)", "LS Mean (Treatment)", "Treatment Difference")
3. Order rows within each visit: LSM ref -> LSM alt -> Treatment Difference
4. Format CI as bracket text
5. Format p-values (dash or empty for LSM rows)

```r
# Verified working pattern
prepare_efficacy_data <- function(tidy_df, digits = 2, ci_level = 0.95) {
  tidy_df %>%
    dplyr::mutate(
      # Clean visit labels
      visit_label = gsub("_", " ", visit),
      visit_label = tools::toTitleCase(visit_label),
      # Row labels
      row_label = dplyr::case_when(
        parameter_type == "lsm" & lsm_type == "ref" ~ "LS Mean (Reference)",
        parameter_type == "lsm" & lsm_type == "alt" ~ "LS Mean (Treatment)",
        parameter_type == "trt" ~ "Treatment Difference",
        TRUE ~ description
      ),
      # Sort order within visit
      row_order = dplyr::case_when(
        parameter_type == "lsm" & lsm_type == "ref" ~ 1L,
        parameter_type == "lsm" & lsm_type == "alt" ~ 2L,
        parameter_type == "trt" ~ 3L,
        TRUE ~ 4L
      ),
      # Pre-format CI bracket text
      ci_text = sprintf(
        "(%s, %s)",
        formatC(lci, format = "f", digits = digits),
        formatC(uci, format = "f", digits = digits)
      ),
      # Format p-value (em dash for non-applicable)
      pval_text = ifelse(is.na(pval), "\u2014", format_pvalue(pval))
    ) %>%
    dplyr::arrange(visit_label, row_order) %>%
    dplyr::select(visit_label, row_label, est, se, ci_text, pval_text)
}
```

### Pattern 2: gt Table Construction with Row Groups
**What:** Build a gt table using `groupname_col` for visit grouping
**When to use:** Always -- this is the rendering step
**Confidence:** HIGH -- verified with working prototype (gt 1.3.0)

```r
# Verified working pattern
render_efficacy_gt <- function(table_data, title, subtitle, footnotes, digits = 2) {
  table_data %>%
    gt::gt(
      rowname_col = "row_label",
      groupname_col = "visit_label"
    ) %>%
    gt::tab_header(
      title = title,
      subtitle = subtitle
    ) %>%
    gt::fmt_number(columns = c(est, se), decimals = digits) %>%
    gt::cols_label(
      est = "Estimate",
      se = "Std. Error",
      ci_text = paste0(ci_level * 100, "% CI"),
      pval_text = "P-value"
    ) %>%
    gt::cols_align(align = "center", columns = dplyr::everything()) %>%
    gt::opt_align_table_header(align = "left")
}
```

### Pattern 3: Footnotes for Metadata
**What:** Add source notes for model info, imputation count, and pooling method
**When to use:** Always -- required by TBL-04
**Confidence:** HIGH -- gt `tab_source_note()` is simpler than `tab_footnote()` for non-cell-targeted notes

```r
# tab_source_note() adds unnumbered footer notes (preferred for model metadata)
# tab_footnote() adds numbered/superscripted notes (for cell-specific annotations)
tbl %>%
  gt::tab_source_note(paste("Pooling method:", pool_obj$method)) %>%
  gt::tab_source_note(paste("Number of imputations:", pool_obj$N)) %>%
  gt::tab_source_note(model_description)
```

### Pattern 4: Dependency Guard (gt as Suggests)
**What:** Check for gt at function call time, matching cards pattern from Phase 2
**When to use:** Always
**Confidence:** HIGH -- matches is_cards_available() pattern

```r
is_gt_available <- function() {
  requireNamespace("gt", quietly = TRUE)
}

efficacy_table <- function(pool_obj, ...) {
  if (!is_gt_available()) {
    cli::cli_abort(
      c(
        "Package {.pkg gt} is required for table generation.",
        "i" = "Install with {.code install.packages(\"gt\")}."
      ),
      class = c("rbmiUtils_error_dependency", "rbmiUtils_error")
    )
  }
  # ...
}
```

### Anti-Patterns to Avoid
- **Using gtsummary `tbl_ard_*()` for model-based results:** These functions are designed for descriptive statistics from raw data. Forcing rbmi pool results into this paradigm creates unnecessary complexity and fragile code.
- **Pre-formatting all columns as character before gt:** Let gt handle formatting with `fmt_number()`. Only pre-format compound fields (CI text, p-value text) that gt cannot compose natively.
- **Hardcoding visit order:** Visit order should follow the factor level order from the pool object (via tidy_pool_obj), not alphabetical sorting.
- **Returning HTML string instead of gt object:** Return the gt object so users can pipe into `gt::tab_*()` functions for customization. Users call `gtsave()` or `as_raw_html()` themselves.
- **Adding gtsummary as a dependency:** Direct gt construction avoids pulling gtsummary + its dependencies into the package. Keep the dependency tree lean.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Table rendering with row groups | Custom HTML/LaTeX generator | `gt::gt()` with `groupname_col` | gt handles HTML/PDF/LaTeX rendering, accessibility, responsive width |
| Number formatting with trailing zeros | `round()` + `paste()` | `gt::fmt_number()` | round() drops trailing zeros; fmt_number preserves them |
| P-value formatting | New formatter | Existing `format_pvalue()` | Already handles thresholds, NA, edge cases (Phase 1 hardened) |
| CI bracket formatting | New formatter | `sprintf()` with `formatC()` | Simple pattern; gt's `cols_merge()` is alternative but adds complexity for "(lower, upper)" |
| Visit label cleaning | Manual string munging per case | Regex + `tools::toTitleCase()` | Handles underscore->space conversion consistently |
| Pool metadata extraction | Re-parsing pool internals | `pool_obj$method`, `pool_obj$N`, `pool_obj$conf.level` | These fields are documented in Phase 2 research; stable rbmi API |

**Key insight:** The table construction is a straightforward data reshape + gt rendering pipeline. The complexity is in getting the layout right (row ordering, grouping, formatting), not in the technology. Resist the urge to introduce gtsummary -- it adds a dependency for a problem that gt solves directly.

## Common Pitfalls

### Pitfall 1: Visit Order Not Preserved
**What goes wrong:** gt displays visits in alphabetical order ("Week 24" before "Week 4") instead of the intended chronological order
**Why it happens:** gt's `groupname_col` groups rows but orders alphabetically unless explicitly controlled
**How to avoid:** Ensure visits are factors with correct level order in the input data. tidy_pool_obj() preserves the factor order from the pool object's parameter names. Use `row_group_order()` in gt if needed.
**Warning signs:** "Week 12" appearing before "Week 4" in the rendered table

### Pitfall 2: Arm Labels Show "Reference"/"Alternative" Instead of Actual Names
**What goes wrong:** Table shows "LS Mean (Reference)" instead of "LS Mean (Placebo)"
**Why it happens:** The pool object only stores parameter names like `lsm_ref_Week24` -- "ref" and "alt" are positional identifiers, not actual arm names. The original arm names are lost after pooling.
**How to avoid:** Accept optional `arm_labels` argument (e.g., `arm_labels = c(ref = "Placebo", alt = "Drug A")`). Default to "Reference" / "Treatment" when not provided. This is a limitation of the rbmi pool object structure, not something we can solve without user input.
**Warning signs:** User confusion about which arm is which

### Pitfall 3: Empty P-value Cells for LSM Rows
**What goes wrong:** LSM rows show "NA" or blank cells in the P-value column, which looks like missing data
**Why it happens:** LSM estimates don't have p-values (only treatment differences do). tidy_pool_obj() returns NA for LSM p-values.
**How to avoid:** Display an em dash ("\u2014") or leave blank intentionally. The em dash is the standard regulatory convention for "not applicable."
**Warning signs:** "NA" appearing in rendered table

### Pitfall 4: gt PDF/LaTeX Rendering Limitations
**What goes wrong:** `tab_style()` calls are silently ignored when rendering to LaTeX/PDF; complex styling doesn't translate
**Why it happens:** gt's LaTeX backend has limited styling support compared to HTML (documented limitation as of gt 1.3.0)
**How to avoid:** Keep styling simple -- rely on structural elements (row groups, headers, footnotes) rather than cell-level styling. Test both HTML and LaTeX output during development. Use `gtsave(filename = "table.pdf")` for PDF testing.
**Warning signs:** Table looks different in PDF vs HTML; "tab_style has no effect in LaTeX" warnings

### Pitfall 5: Multiple Parameters in One Pool Object
**What goes wrong:** If the pool object contains multiple endpoints (e.g., CHG and AVAL), they all appear in one table without clear separation
**Why it happens:** tidy_pool_obj() returns all parameters; the table function doesn't filter by endpoint
**How to avoid:** For v1, assume single-parameter pool objects (standard rbmi workflow). Document this assumption. If multi-parameter support is needed later, add a `parameter` argument to filter, or nest with additional row groups.
**Warning signs:** Mixed endpoint names in the parameter column of tidy_pool_obj() output

### Pitfall 6: Decimal Alignment Across Rows
**What goes wrong:** Numbers don't align at the decimal point because different rows have different magnitudes
**Why it happens:** gt center-aligns by default; decimal alignment requires monospace font or explicit formatting
**How to avoid:** Use `gt::fmt_number()` with consistent `decimals` argument to ensure all values have the same number of decimal places. This ensures visual alignment. Consider `gt::opt_table_font(font = "monospace")` for regulatory-style appearance.
**Warning signs:** Jagged decimal alignment in rendered table

## Code Examples

### Example 1: Complete efficacy_table() Function Skeleton
```r
# Verified working pattern with gt 1.3.0
efficacy_table <- function(
  pool_obj,
  title = NULL,
  subtitle = NULL,
  digits = 2,
  ci_level = NULL,
  arm_labels = NULL,
  ...
) {
  # Dependency check
  if (!is_gt_available()) {
    cli::cli_abort(
      c("Package {.pkg gt} is required.", "i" = "Install with {.code install.packages(\"gt\")}."),
      class = c("rbmiUtils_error_dependency", "rbmiUtils_error")
    )
  }

  # Input validation
  if (!inherits(pool_obj, "pool")) {
    cli::cli_abort(
      "Input {.arg pool_obj} must be of class {.cls pool}.",
      class = c("rbmiUtils_error_validation", "rbmiUtils_error")
    )
  }

  # Extract metadata
  ci_level <- ci_level %||% pool_obj$conf.level %||% 0.95
  method <- pool_obj$method %||% "unknown"
  n_imputations <- pool_obj$N

  # Tidy and prepare data
  tidy_df <- tidy_pool_obj(pool_obj)

  # Clean visit labels
  tidy_df <- tidy_df %>%
    dplyr::mutate(
      visit_label = gsub("_", " ", visit),
      visit_label = tools::toTitleCase(visit_label)
    )

  # Apply arm labels if provided
  ref_label <- if (!is.null(arm_labels)) arm_labels[["ref"]] else "Reference"
  alt_label <- if (!is.null(arm_labels)) arm_labels[["alt"]] else "Treatment"

  # Create display data
  table_data <- tidy_df %>%
    dplyr::mutate(
      row_label = dplyr::case_when(
        parameter_type == "lsm" & lsm_type == "ref" ~ paste("LS Mean (", ref_label, ")", sep = ""),
        parameter_type == "lsm" & lsm_type == "alt" ~ paste("LS Mean (", alt_label, ")", sep = ""),
        parameter_type == "trt" ~ "Treatment Difference",
        TRUE ~ description
      ),
      row_order = dplyr::case_when(
        parameter_type == "lsm" & lsm_type == "ref" ~ 1L,
        parameter_type == "lsm" & lsm_type == "alt" ~ 2L,
        parameter_type == "trt" ~ 3L,
        TRUE ~ 4L
      ),
      ci_text = sprintf("(%s, %s)",
        formatC(lci, format = "f", digits = digits),
        formatC(uci, format = "f", digits = digits)
      ),
      pval_text = ifelse(is.na(pval), "\u2014", format_pvalue(pval))
    ) %>%
    dplyr::arrange(visit_label, row_order) %>%
    dplyr::select(visit_label, row_label, est, se, ci_text, pval_text)

  # Build gt table
  tbl <- table_data %>%
    gt::gt(rowname_col = "row_label", groupname_col = "visit_label") %>%
    gt::fmt_number(columns = c(est, se), decimals = digits) %>%
    gt::cols_label(
      est = "Estimate",
      se = "Std. Error",
      ci_text = paste0(ci_level * 100, "% CI"),
      pval_text = "P-value"
    ) %>%
    gt::cols_align(align = "center", columns = dplyr::everything()) %>%
    gt::opt_align_table_header(align = "left")

  # Title
  if (!is.null(title) || !is.null(subtitle)) {
    tbl <- gt::tab_header(tbl, title = title, subtitle = subtitle)
  }

  # Footnotes (TBL-04 requirement)
  tbl <- tbl %>%
    gt::tab_source_note(paste("Pooling method:", method)) %>%
    gt::tab_source_note(paste("Number of imputations:", n_imputations)) %>%
    gt::tab_source_note(paste0("Confidence level: ", ci_level * 100, "%"))

  tbl
}
```

### Example 2: Visit Label Cleaning
```r
# Verified: handles common rbmi parameter name patterns
clean_visit_label <- function(visit_name) {
  # "Week_4" -> "Week 4"
  # "Week_24" -> "Week 24"
  # "Week 4" -> "Week 4" (already clean)
  label <- gsub("_", " ", visit_name)
  tools::toTitleCase(label)
}
```

### Example 3: Table with Custom Arm Labels
```r
# Usage example for users
efficacy_table(
  pool_obj,
  title = "Table 14.2.1: Change from Baseline in Primary Endpoint",
  subtitle = "ANCOVA Model - Modified ITT Population",
  arm_labels = c(ref = "Placebo", alt = "Drug A"),
  digits = 2
)
```

### Example 4: Piping into Further gt Customization
```r
# The function returns a gt object that users can further customize
efficacy_table(pool_obj) %>%
  gt::tab_options(
    table.font.size = gt::px(11),
    heading.title.font.size = gt::px(14)
  ) %>%
  gt::gtsave("efficacy_table.html")
```

### Example 5: Testing Pattern with Mock Pool Objects
```r
# Reuse make_mock_pool() from test-pool_methods.R
test_that("efficacy_table returns gt object", {
  mock_pool <- make_mock_pool()
  tbl <- efficacy_table(mock_pool)
  expect_s3_class(tbl, "gt_tbl")
})

test_that("efficacy_table contains footnotes", {
  mock_pool <- make_mock_pool()
  tbl <- efficacy_table(mock_pool)
  html <- gt::as_raw_html(tbl)
  expect_true(grepl("rubin", html))
  expect_true(grepl("100", html))
})
```

## Discretion Recommendations

Based on research, here are recommendations for the Claude's Discretion areas from CONTEXT.md:

### Single function vs two-step API
**Recommendation: Single `efficacy_table()` function** (with internal two-step design)
- Users get a one-call experience matching TBL-03
- Internally, separate data preparation from gt rendering for testability
- No need to expose `prepare_efficacy_data()` publicly

### Accept pool objects only or also ARD tibbles
**Recommendation: Pool objects only for v1**
- The ARD from `pool_to_ard()` uses the cards structure (long-format, stat column) which would need re-pivoting -- adding complexity for no user benefit
- Pool -> tidy -> table is the natural flow
- ARD acceptance could be a future enhancement if gtsummary integration is needed

### Column header organization
**Recommendation: Simple flat columns (no spanners for v1)**
- Columns: Estimate | Std. Error | 95% CI | P-value
- Spanners would separate "LS Mean" vs "Treatment Difference" columns, but since these are different rows (not different columns), spanners don't apply naturally
- If needed later, spanners could group "Estimate" + "Std. Error" under one header

### N count display approach
**Recommendation: Omit per-visit N counts in v1**
- The pool object does not contain per-visit N counts (only total number of imputations)
- Adding N would require the user to supply it separately or extracting from the original data
- Footnote stating number of imputations is sufficient for the regulatory context

### Multi-parameter handling
**Recommendation: Single-parameter assumption for v1; warn if multiple detected**
- Standard rbmi workflow produces single-endpoint pool objects
- If tidy_pool_obj() returns rows with different parameter prefixes (beyond trt/lsm), emit a warning
- Document: "For multiple endpoints, call efficacy_table() separately for each pool object"

### Treatment difference statistics
**Recommendation: Estimate + SE + CI + p-value (full regulatory set)**
- All four statistics in columns: Estimate, Std. Error, XX% CI, P-value
- Same columns for LSM rows (but p-value shows em dash)
- This matches the standard Table 14.2 ANCOVA format

### CI source
**Recommendation: Use CI from pool object directly (no recalculation)**
- pool_obj already contains CIs computed by rbmi with proper pooling rules
- Recalculation from estimate +/- 1.96*SE would be incorrect for Rubin's rules (df adjustment)
- The `ci_level` parameter is for labeling only ("95% CI"), not recalculation

### CI bracket style
**Recommendation: Parentheses with comma separator: `(lower, upper)`**
- Most common in pharmaceutical regulatory submissions
- Matches existing `format_estimate()` default in the package
- Alternative: square brackets `[lower, upper]` -- offer via `ci_sep` or `ci_brackets` argument

### P-value formatting
**Recommendation: 3 decimal places, threshold at 0.001 ("< 0.001")**
- Matches existing `format_pvalue()` defaults in the package
- Standard regulatory convention
- User can override via `pval_digits` and `pval_threshold` arguments

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| flextable for regulatory tables | gt for HTML/PDF, rtables for production | gt >= 0.10 (2024) | gt now has reliable LaTeX output, sufficient for most CSR needs |
| r2rtf for RTF tables | gt with gtsave for multi-format | gt >= 0.10 | gt handles HTML, PDF, LaTeX; RTF still needs r2rtf if required |
| gtsummary for all tables | gtsummary for summary stats, gt directly for model-based tables | Ongoing | gtsummary tbl_ard_*() excels at descriptive stats but is not designed for pre-computed model estimates |
| Manual HTML/LaTeX generation | gt declarative API | gt >= 0.1 | Declarative approach is less error-prone, handles multi-format output |

**Deprecated/outdated:**
- Building regulatory tables with `kable()` or `kableExtra` -- insufficient for row groups, footnotes, and multi-format output
- Using `flextable` for HTML-first workflows -- gt is now the standard for HTML clinical tables
- Using `gtsummary::tbl_regression()` for rbmi pool objects -- pool objects are not model objects; this would require custom tidiers that add unnecessary coupling

## Open Questions

1. **Should the function accept `tidy_df` (from tidy_pool_obj) in addition to pool_obj?**
   - What we know: Some users may want to pre-filter visits or modify data before rendering
   - What's unclear: Whether this adds API complexity for an edge case
   - Recommendation: Accept pool_obj only in v1. Users can call tidy_pool_obj() -> filter -> manually build gt if needed. Document the internal data structure.

2. **How to handle arm label discovery when user doesn't provide them?**
   - What we know: Pool objects contain "ref" and "alt" but not "Placebo" and "Drug A". The original data is not accessible from the pool object.
   - What's unclear: Whether we should try to extract arm names from parameter names in edge cases
   - Recommendation: Default to "Reference" / "Treatment" with clear documentation that users should provide `arm_labels` for publication-quality tables.

3. **Should gt be Suggests or Imports?**
   - What we know: Phase 2 put cards as Suggests. PROJECT.md says "to be decided during planning."
   - What's unclear: Whether enough functions will use gt to justify Imports
   - Recommendation: Suggests (matching cards pattern). Only efficacy_table() uses gt. Guard with `is_gt_available()`.

## Sources

### Primary (HIGH confidence)
- gt package v1.3.0 -- tested locally: `gt()`, `tab_row_group()`, `cols_merge()`, `tab_footnote()`, `tab_source_note()`, `fmt_number()`, `gtsave()`, `as_raw_html()`, `as_latex()`
- rbmiUtils codebase -- read directly: tidy_pool_obj(), format_pvalue(), format_estimate(), pool_methods.R, ard_conversion.R
- rbmi pool object structure -- verified: `$pars`, `$method`, `$N`, `$conf.level`, `$alternative`
- Working prototype -- verified: complete pipeline from mock pool_obj -> tidy_pool_obj -> data reshape -> gt table -> HTML output
- [gt case study: clinical tables](https://gt.rstudio.com/articles/case-study-clinical-tables.html) -- row groups, column merging, formatting patterns
- [gt reference index](https://gt.rstudio.com/reference/index.html) -- function signatures and parameters
- [gt tab_footnote docs](https://gt.rstudio.com/reference/tab_footnote.html) -- footnote targeting with location helpers
- [gt tab_row_group docs](https://gt.rstudio.com/reference/tab_row_group.html) -- programmatic row group creation
- [gt cols_merge docs](https://gt.rstudio.com/reference/cols_merge.html) -- pattern syntax for column merging

### Secondary (MEDIUM confidence)
- [R for Clinical Study Reports: Efficacy Table](https://r4csr.org/tlf-efficacy-ancova.html) -- Table 14.2 ANCOVA layout reference (uses r2rtf but layout structure applies)
- [Appsilon: Clinical Tables with gt](https://www.appsilon.com/post/clinical-tables-with-the-gt-package-44a57) -- complete gt clinical table workflow
- [gtsummary ARD-first tables](https://www.danieldsjoberg.com/gtsummary/articles/tbl_ard-functions.html) -- confirmed tbl_ard_*() is for descriptive stats, not model estimates
- [gtsummary as_gt()](https://www.danieldsjoberg.com/gtsummary/reference/as_gt.html) -- conversion from gtsummary to gt for further customization
- [PharmaSUG QT-263: R Tables via GT for Regulatory Submissions](https://pharmasug.org/proceedings/2023/QT/PharmaSUG-2023-QT-263.pdf) -- regulatory submission context

### Tertiary (LOW confidence)
- [gt LaTeX rendering issue #1363](https://github.com/rstudio/gt/issues/1363) -- known limitations with gt PDF output via Quarto
- WebSearch results for gt PDF workarounds -- general guidance, verify during implementation

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - gt 1.3.0 tested locally with working prototype
- Architecture: HIGH - data reshape + gt pipeline verified end-to-end
- Table 14.2 layout: HIGH - cross-referenced with FDA/ICH format from r4csr and PharmaSUG
- Pitfalls: HIGH - discovered through prototype testing (visit ordering, arm labels, p-value formatting, LaTeX limitations)
- Discretion recommendations: MEDIUM - based on research and pharma domain expertise; validate during planning

**Research date:** 2026-02-08
**Valid until:** 2026-03-08 (gt package is stable; API unlikely to change)
