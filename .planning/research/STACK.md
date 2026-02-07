# Technology Stack: Clinical Trial Reporting Layer for rbmiUtils

**Project:** rbmiUtils reporting milestone (ARD/gtsummary integration)
**Researched:** 2026-02-07
**Overall confidence:** MEDIUM-HIGH (versions verified via CRAN; integration patterns verified via official docs and workshop materials)

---

## Recommended Stack

### Core ARD Layer

| Technology | Version | Purpose | Why | Confidence |
|------------|---------|---------|-----|------------|
| cards | 0.7.1 (2025-12-02) | ARD data structure, `ard_identity()`, `as_card()`, `bind_ard()` | The canonical CDISC ARD implementation. Co-developed by Roche, GSK, Novartis. 52k+ downloads/month. Provides the standard data structure that gtsummary consumes. | HIGH |
| cardx | 0.3.2 (2026-02-05) | Extended ARD functions: `ard_emmeans_emmeans()`, `ard_emmeans_contrast()`, `ard_stats_t_test()` | Provides LS mean and contrast ARD functions via emmeans integration. Directly relevant for rbmi pool objects that contain LS means and treatment comparisons. Soft-depends on emmeans (not installed by default). | HIGH |
| gtsummary | 2.5.0 (2025-12-05) | Table generation from ARDs via `tbl_ard_*()` family | The dominant R package for clinical summary tables. ARD-first workflow (`tbl_ard_summary`, `tbl_ard_wide_summary`, `tbl_ard_continuous`) is the endorsed approach for 2025+. Imports both cards and cardx. | HIGH |
| gt | 1.3.0 (2026-01-22) | Table rendering to HTML, RTF, LaTeX/PDF, Word | Default rendering backend for gtsummary. Supports `as_rtf()`, `as_latex()`, `as_word()`, `gtsave()`. Active Posit development with pharma-specific features. | HIGH |

### Visualization Layer

| Technology | Version | Purpose | Why | Confidence |
|------------|---------|---------|-----|------------|
| ggplot2 | >= 3.5.0 | Base plotting engine for forest plots, responder bar charts | Already in rbmiUtils Suggests. Universal foundation; all clinical viz packages build on it. | HIGH |
| forestploter | 1.1.3 (2025-04-13) | Forest plots for LS mean treatment differences with CIs | Treats forest plot as a table (rows/columns alignment). Full control over layout. Supports multiple CI columns and grouped display. Best fit for rbmi-style visit-by-visit treatment comparisons. | MEDIUM |

### Supporting Libraries

| Library | Version | Purpose | When to Use | Confidence |
|---------|---------|---------|-------------|------------|
| gtreg | 0.4.2 (2025-11-25) | Regulatory table templates (AE tables, reg summaries) | Only if adverse event tables or `tbl_reg_summary()` are needed. Built on gtsummary. | MEDIUM |
| flextable | >= 6.0 | Alternative rendering backend for Word/PowerPoint | When `as_gt()` output does not meet specific Word formatting requirements. gtsummary supports `as_flex_table()`. | MEDIUM |
| r2rtf | >= 1.1 | Low-level RTF generation | Only if gt's `as_rtf()` does not meet pixel-perfect RTF requirements for a specific submission. More control but more manual work. | LOW |

---

## The ARD Data Structure

The ARD (Analysis Results Dataset) is the central data contract between analysis and reporting. Understanding its column structure is essential for the rbmiUtils bridge.

### Required Columns

| Column | Type | Description |
|--------|------|-------------|
| `variable` | character | Name of the analyzed variable (e.g., "CHG") |
| `stat_name` | character | Statistic identifier (e.g., "estimate", "std.error", "p.value") |
| `stat_label` | character | Human-readable label (e.g., "Estimate", "Standard Error", "P-value") |
| `stat` | list | The calculated value (stored as list column to support mixed types) |

### Grouping Columns (Optional, used for stratified results)

| Column | Type | Description |
|--------|------|-------------|
| `group1` | character | First grouping variable name (e.g., "TRT") |
| `group1_level` | list | Level of group1 (e.g., "Drug A") |
| `group2` | character | Second grouping variable name (e.g., "AVISIT") |
| `group2_level` | list | Level of group2 (e.g., "Week 24") |
| `variable_level` | list | Level of the variable (for categorical) |

### Metadata Columns

| Column | Type | Description |
|--------|------|-------------|
| `context` | character | Analysis context label (e.g., "identity", "emmeans") |
| `fmt_fn` | list | Formatting function for display |
| `warning` | list | Warnings captured during calculation |
| `error` | list | Errors captured during calculation |

### Key Insight for rbmiUtils

The `stat` column is a **list column** -- each element can be numeric, character, or any R object. This is how ARDs store heterogeneous results (estimates, p-values, labels) in a single tidy data frame. The `ard_identity()` function is the primary entry point for wrapping pre-calculated statistics into this structure.

---

## How rbmi Pool Objects Map to ARDs

### Current tidy_pool_obj() Output

```
parameter | description | visit | parameter_type | lsm_type | est | se | lci | uci | pval
```

### Mapping Strategy: Use `ard_identity()` from cards

The `ard_identity()` function accepts pre-calculated statistics (as a named list) and reformats them into the ARD structure. This is the correct approach for rbmi pool objects because:

1. **rbmi already computed the statistics** -- we are not re-running analysis, just reformatting results
2. **The pool object has a fixed schema** -- est, se, lci, uci, pval are already computed
3. **ard_identity() preserves values exactly** -- no recomputation, no rounding

### Proposed Conversion Pattern

For each row of `tidy_pool_obj()` output, create an ARD with:

| tidy_pool_obj column | ARD column | Mapping |
|---------------------|------------|---------|
| `parameter_type` + `lsm_type` | `group1` / `group1_level` | Treatment arm or comparison context |
| `visit` | `group2` / `group2_level` | Visit label |
| `"CHG"` (or configurable) | `variable` | The outcome variable |
| `est` | `stat` where `stat_name = "estimate"` | Point estimate |
| `se` | `stat` where `stat_name = "std.error"` | Standard error |
| `lci` | `stat` where `stat_name = "conf.low"` | Lower CI bound |
| `uci` | `stat` where `stat_name = "conf.high"` | Upper CI bound |
| `pval` | `stat` where `stat_name = "p.value"` | P-value |

Each row in tidy_pool_obj produces 5 rows in the ARD (one per statistic). The `bind_ard()` function from cards combines them.

### Alternative: `as_card()` for Direct Conversion

If the tidy_pool_obj output can be reshaped to have `stat_name`, `stat_label`, and `stat` columns, `as_card()` converts the data frame directly to an ARD of class `card`. This may be simpler for batch conversion.

### Why NOT cardx's ard_emmeans_*()

The `ard_emmeans_emmeans()` and `ard_emmeans_contrast()` functions in cardx operate on regression model objects (lm, glm, etc.) and call emmeans internally. They are designed for computing LS means from scratch. rbmi's pool objects are **already pooled results** from multiple imputations -- you cannot pass them to emmeans. Therefore, `ard_identity()` (wrapping pre-computed results) is the correct approach, not `ard_emmeans_*()` (recomputing from a model).

---

## gtsummary ARD-First Table Functions

### Function Selection Guide

| Function | Use Case | Relevant for rbmiUtils? |
|----------|----------|------------------------|
| `tbl_ard_summary()` | Descriptive stats (N, mean, median) by group | YES -- for baseline characteristics if added |
| `tbl_ard_wide_summary()` | Stats in columns (e.g., LS Mean columns per arm) | YES -- primary function for efficacy results tables |
| `tbl_ard_continuous()` | Continuous variable summary stratified by group | MAYBE -- for continuous endpoint summaries |
| `tbl_ard_hierarchical()` | Nested/hierarchical data | NO -- not relevant for rbmi results |
| `tbl_ard_strata()` | Stratified tables from ARDs | MAYBE -- for subgroup analyses |

### Recommended Primary Function: `tbl_ard_wide_summary()`

For rbmi efficacy results (LS means and treatment comparisons), `tbl_ard_wide_summary()` is the best fit because:

1. Clinical efficacy tables typically show statistics in columns (Estimate, CI, P-value)
2. Treatment arms appear as column groups
3. Visits appear as rows
4. This matches the standard regulatory table layout

### Workflow

```
rbmi::pool() --> tidy_pool_obj() --> as_ard() [new function] --> tbl_ard_wide_summary() --> as_gt() --> gtsave()
```

---

## gt Rendering Capabilities for Regulatory Submissions

### Output Format Support

| Format | Function | Regulatory Use | Quality |
|--------|----------|---------------|---------|
| HTML | `as_raw_html()` / `gtsave(".html")` | Internal review, interactive reports | Excellent |
| RTF | `as_rtf()` / `gtsave(".rtf")` | FDA submission appendices | Good -- may need `pharmaRTF` wrapper for pixel-perfect control |
| LaTeX/PDF | `as_latex()` / `gtsave(".pdf")` | CSR (Clinical Study Report) appendices | Good -- requires longtable, booktabs LaTeX packages |
| Word (.docx) | `as_word()` / `gtsave(".docx")` | Medical writing, reviewer comments | Good |

### Key gt Features for Clinical Tables

- **`tab_header()`**: Title and subtitle (required for regulatory tables)
- **`tab_footnote()`**: Footnotes with locator functions
- **`tab_source_note()`**: Data source annotation
- **`tab_spanner()`**: Spanning column headers (treatment arms)
- **`fmt_number()`**: Decimal precision control
- **`cols_align()`**: Column alignment (typically decimal-aligned)
- **`tab_style()`**: Conditional formatting

### Recommendation

Use gt as the primary rendering engine. Fall back to flextable (`as_flex_table()`) only for Word documents requiring specific formatting that gt cannot produce. Use r2rtf only as a last resort for submission-critical RTF formatting.

---

## Visualization Stack

### Forest Plots

**Recommendation: forestploter (1.1.3)**

| Feature | forestploter | ggforestplot | forester | Custom ggplot2 |
|---------|-------------|--------------|----------|----------------|
| Table-aligned layout | YES (native) | NO | YES | Manual |
| Multiple CI columns | YES | NO | NO | Manual |
| Grouped by visit | YES | NO | YES | Manual |
| Publication quality | HIGH | MEDIUM | HIGH | Varies |
| CRAN maintained | YES (2025) | Last update 2022 | Not on CRAN | N/A |
| Customization | HIGH (theme system) | LOW | MEDIUM | Unlimited |

**Rationale:** forestploter treats the forest plot as a table where CIs are rendered in specific columns. This is the standard layout for regulatory forest plots showing treatment differences across visits. It is actively maintained on CRAN (April 2025). ggforestplot has not been updated since 2022.

**Alternative for interactive use:** forestly (0.1.3, 2025-09-10) provides interactive forest plots using plotly. Suitable for Shiny apps or interactive HTML reports but not for static regulatory documents.

### Responder Bar Charts

**Recommendation: ggplot2 with `geom_col()` + `geom_errorbar()`**

No specialized package needed. Standard ggplot2 bar charts with:
- `geom_col()` for proportion bars per treatment arm
- `geom_errorbar()` for confidence intervals
- `facet_wrap()` for visit-by-visit display
- `scale_y_continuous(labels = scales::percent)` for percentage axis

This is simpler and more maintainable than adding a dependency for what amounts to a basic bar chart.

---

## What NOT to Use (and Why)

| Package | Why Not | Use Instead |
|---------|---------|-------------|
| rtables | Different ecosystem (NEST/Roche). Incompatible paradigm with cards/gtsummary. Steep learning curve. | gtsummary + cards |
| tfrmt | Requires a separate metadata layer. Adds complexity without clear benefit when gtsummary already consumes ARDs. | gtsummary `tbl_ard_*()` |
| chevron | Built on rtables. Same ecosystem concern. | gtsummary |
| Tplyr | Older approach. Does not integrate with ARD standard. | cards + gtsummary |
| tidytlg | SAS-migration tool. Not designed for ARD-first workflows. | cards + gtsummary |
| huxtable | Generalist table package. No clinical-specific features. gt is more capable and better maintained. | gt |
| ggforestplot | Last updated 2022. Vertical layout only. Not maintained. | forestploter |
| visR | Focused on survival analysis (KM plots). Not relevant for MMRM/ANCOVA results. Last CRAN update 2024-03. | ggplot2 + forestploter |

---

## Dependency Management Strategy

### New Imports (required)

```
cards (>= 0.7.0)
```

**Rationale:** Only cards is needed as a hard dependency because `ard_identity()`, `as_card()`, and `bind_ard()` are the core bridge functions. The cards package has minimal dependencies itself (dplyr, tidyr, rlang, cli, glue, lifecycle, tidyselect) -- all of which rbmiUtils already imports or depends on transitively.

### New Suggests (optional, loaded on use)

```
gtsummary (>= 2.0.0),
gt (>= 1.0.0),
cardx (>= 0.2.0),
forestploter (>= 1.0.0)
```

**Rationale:** Keep these as Suggests because:
1. Not all users need table rendering (some only need the ARD for downstream use)
2. gtsummary pulls in gt, cardx, and many other packages -- heavy install
3. forestploter is only needed for forest plot generation
4. Following the same pattern as cardx itself (soft dependencies on statistical packages)

### Installation

```r
# Core (minimal)
install.packages("cards")

# Full reporting stack
install.packages(c("gtsummary", "gt", "forestploter"))

# cardx comes automatically with gtsummary
```

---

## Version Compatibility Matrix

| Package | Min Version | Reason |
|---------|-------------|--------|
| R | >= 4.2 | Required by cardx and gtsummary |
| cards | >= 0.7.0 | `ard_identity()` API stable from this version |
| cardx | >= 0.2.0 | emmeans ARD functions available |
| gtsummary | >= 2.0.0 | ARD-first `tbl_ard_*()` functions introduced in v2.x rewrite |
| gt | >= 1.0.0 | Stable `as_rtf()`, `as_word()` APIs |
| dplyr | >= 1.1.0 | Already required by rbmiUtils |

**IMPORTANT:** gtsummary 2.x was a major rewrite that introduced the ARD-first workflow. Versions 1.x do NOT have `tbl_ard_*()` functions. The minimum version must be >= 2.0.0.

---

## Confidence Assessment

| Recommendation | Confidence | Basis |
|---------------|------------|-------|
| cards as ARD foundation | HIGH | CRAN verified (v0.7.1, Dec 2025). Official CDISC ARD implementation. Backed by Roche/GSK/Novartis. |
| ard_identity() for pool conversion | HIGH | Official docs confirm it wraps pre-calculated stats into ARD format. Exact use case. |
| gtsummary tbl_ard_wide_summary() | MEDIUM | Function exists and is documented, but marked "experimental". May change API. |
| gt for multi-format rendering | HIGH | CRAN verified (v1.3.0, Jan 2026). as_rtf/as_word/as_latex all documented and stable. |
| forestploter for forest plots | MEDIUM | CRAN maintained (Apr 2025). Good fit for clinical forest plots. Less community adoption data available. |
| ggplot2 for responder charts | HIGH | Standard approach. No specialized package needed. |
| as_card() for batch conversion | MEDIUM | Function exists but limited examples in docs. May need experimentation. |

---

## Sources

### CRAN Package Pages (version verification)
- [cards 0.7.1 on CRAN](https://cran.r-project.org/web/packages/cards/index.html)
- [cardx 0.3.2 on CRAN](https://cran.r-project.org/web/packages/cardx/index.html)
- [gtsummary 2.5.0 on CRAN](https://cran.r-project.org/web/packages/gtsummary/index.html)
- [gt 1.3.0 on CRAN](https://cran.r-project.org/web/packages/gt/index.html)
- [forestploter 1.1.3 on CRAN](https://cran.r-project.org/web/packages/forestploter/index.html)
- [gtreg 0.4.2 on CRAN](https://cran.r-project.org/web/packages/gtreg/index.html)
- [forestly 0.1.3 on CRAN](https://cran.r-project.org/web/packages/forestly/index.html)

### Official Documentation
- [cards package docs](https://insightsengineering.github.io/cards/main/index.html)
- [cardx package docs](https://insightsengineering.github.io/cardx/main/)
- [gtsummary ARD-first tables](https://www.danieldsjoberg.com/gtsummary/articles/tbl_ard-functions.html)
- [gtsummary tbl_ard_wide_summary](https://www.danieldsjoberg.com/gtsummary/reference/tbl_ard_wide_summary.html)
- [gt as_rtf reference](https://gt.rstudio.com/reference/as_rtf.html)
- [cards ard_identity reference](https://insightsengineering.github.io/cards/main/reference/ard_identity.html)
- [cards as_card reference](https://insightsengineering.github.io/cards/main/reference/as_card.html)
- [forestploter vignette](https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html)

### Workshop/Presentation Materials
- [CDISC COSA Spotlight: ARD-based Reporting with cards and gtsummary (2025)](https://www.danieldsjoberg.com/CDISC-COSA-Spotlight-ARD-gtsummary-2025/)
- [ARD PHUSE Workshop 2025: cards overview](https://www.danieldsjoberg.com/ARD-PHUSE-workshop-2025/slides/01-intro-cards/)
- [ARD PHUSE Workshop 2025: cardx extras](https://www.danieldsjoberg.com/ARD-PHUSE-workshop-2025/slides/02-intro-cardx/)
- [Posit Conf 2025: End-to-End Pharmaverse with gtsummary](https://posit-conf-2025.github.io/pharmaverse/slides/06-tables-gtsummary/)
- [R Consortium: Supercharging Statistical Analysis with ARDs](https://r-consortium.org/posts/supercharging-statistical-analysis-with-ards-and-the-cards-r-package/)

### Ecosystem References
- [pharmaverse.org TLG packages](https://pharmaverse.org/e2eclinical/tlg/)
- [R Consortium RTRS Working Group](https://rconsortium.github.io/rtrs-wg/)
- [R for Clinical Study Reports (r4csr.org)](https://r4csr.org/tlf-overview.html)
