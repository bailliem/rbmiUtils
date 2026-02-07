# Feature Landscape

**Domain:** Clinical trial reporting utilities for rbmi (reference-based multiple imputation)
**Researched:** 2026-02-07
**Overall confidence:** MEDIUM-HIGH (based on codebase inspection + ecosystem documentation + pharmaverse standards)

## Context: What rbmiUtils Has Today

The package currently provides a pipeline from pooled MI results to formatted strings:

```
rbmi::pool() -> tidy_pool_obj() -> tibble(est, se, lci, uci, pval)
                                      |
                                      +-> format_results() -> formatted strings
                                      +-> extract_trt_effects() / extract_lsm() -> filtered rows
                                      +-> combine_results() -> merged analyses
```

The gap: These formatted tibbles are "almost there" but cannot flow into the pharmaverse ARD/table ecosystem, produce regulatory-quality tables, or generate standard clinical trial figures. Users must manually wrangle from `format_results()` output into their table framework of choice.

---

## Table Stakes

Features users expect from an rbmi extension that claims to produce regulatory tables and figures. Missing any of these makes the package feel incomplete for its stated purpose.

| # | Feature | Why Expected | Complexity | Dependencies | Notes |
|---|---------|-------------|------------|--------------|-------|
| TS-1 | **ARD conversion from tidy pool results** | The pharmaverse has standardized on ARD (Analysis Results Datasets) via the {cards} package as the interchange format between analysis and reporting. Any package producing clinical trial results needs ARD output to be interoperable. | Medium | cards (>= 0.3.0) | Use `cards::tidy_as_ard()` or build a custom `as_ard.pool()` / `pool_to_ard()` function. The tidy pool tibble columns (est, se, lci, uci, pval, parameter, visit, parameter_type) map naturally to ARD's group1/variable/stat_name/stat structure. |
| TS-2 | **Efficacy summary table (regulatory style)** | Every clinical study report (ICH E3 Section 11) includes a primary efficacy table showing LS means by arm, treatment differences, CIs, and p-values by visit. This is the most requested output from any MI analysis package. | Medium-High | gt or gtsummary, cards | Standard layout: rows = visits, columns = Reference LS Mean, Treatment LS Mean, Difference (95% CI), P-value. Must support both single-endpoint and multi-endpoint layouts. |
| TS-3 | **Forest plot for treatment effects** | Standard visualization in clinical trial reports showing treatment effects across visits or subgroups with point estimates and CIs. Expected in every regulatory submission with longitudinal efficacy data. | Medium | ggplot2 | Build with ggplot2 (already Suggests dependency). Input: tidy pool results from `extract_trt_effects()`. Output: ggplot object. Must handle visit-by-visit and subgroup layouts. |
| TS-4 | **Print/summary methods for pool objects** | rbmi's `print.pool()` shows raw numeric output with minimal formatting. Users expect enhanced print methods that show rounded estimates, formatted CIs, and clear parameter labeling -- the equivalent of what `tidy_pool_obj()` produces but at the console. | Low | None (base R S3) | This is an S3 method wrapping existing tidy/format logic. Can leverage `tidy_pool_obj()` + `format_results()` internally. Register as `print.pool` or provide a `print_pool()` wrapper to avoid masking rbmi's method. |
| TS-5 | **Responder analysis bar chart** | Binary endpoint reporting requires grouped bar charts showing proportion responding by treatment arm across visits. This is the visual counterpart to `gcomp_responder_multi()` results. | Low-Medium | ggplot2 | Input: pooled responder results (from gcomp pathway). Output: grouped bar chart with error bars (CIs). Standard pharma reporting visualization. |
| TS-6 | **Column formatting controls for tables** | Regulatory tables require specific decimal precision, CI bracket styles, p-value thresholds, and column ordering. The existing `format_pvalue()`, `format_estimate()`, and `format_results_table()` cover basics but need to integrate with gt/gtsummary theming. | Low | gt or gtsummary | Extend existing formatters to produce gt-compatible output. May need a `theme_rbmi()` or formatting configuration object. |

### Dependency Note on TS-1 and TS-2

TS-1 (ARD) is a prerequisite for TS-2 (efficacy table) if the table is built via the gtsummary `tbl_ard_*()` pathway. However, TS-2 can alternatively be built directly with gt without ARD. **Recommendation:** Build TS-1 first, then TS-2 consumes the ARD.

---

## Differentiators

Features that set rbmiUtils apart from generic clinical trial reporting packages. Not expected, but create significant value for rbmi users specifically.

| # | Feature | Value Proposition | Complexity | Dependencies | Notes |
|---|---------|-------------------|------------|--------------|-------|
| D-1 | **One-call efficacy table from pool object** | `efficacy_table(pool_obj)` that goes directly from `rbmi::pool()` output to a gt table with no intermediate steps. Competing workflows (tern/rtables, gtsummary) require 5-10 lines of configuration. A single function that "just works" for the rbmi use case is a strong differentiator. | Medium | gt, cards | Internally chains: pool -> tidy -> ARD -> gt table. The value is in the opinionated defaults: visit rows, arm columns, standard regulatory formatting. Power users can still use the intermediate steps. |
| D-2 | **Sensitivity analysis comparison table** | Side-by-side table comparing primary analysis vs. tipping point / delta-adjusted results. This is unique to the MI workflow -- rbmi users always run sensitivity analyses and need to present them together. No existing package does this for rbmi output. | Medium | gt | Leverages `combine_results()` to stack analyses, then generates a comparison layout. Columns: Primary Est (95% CI), Sensitivity Est (95% CI), Difference. |
| D-3 | **Integrated ARD with MI metadata** | Standard ARD captures point estimates but not MI-specific metadata (number of imputations, pooling method, fraction of missing information). An rbmi-aware ARD function can attach this metadata, making the ARD self-documenting for regulatory review. | Low-Medium | cards | Extend ARD with additional stat_name entries for n_imputations, pooling_method, fmi (fraction of missing information) if available from pool object. |
| D-4 | **Forest plot with sensitivity overlay** | Forest plot showing both primary and sensitivity analysis estimates side-by-side for each visit/subgroup. Unique to MI workflows where sensitivity analyses are mandatory. | Medium | ggplot2 | Extend TS-3 to accept multiple analysis results (e.g., from `combine_results()`). Display as offset points or dual-panel layout. |
| D-5 | **Responder analysis with treatment difference annotations** | Bar chart from TS-5 with annotated treatment differences, CIs, and p-values directly on the plot. Makes the figure self-contained for regulatory documents. | Low-Medium | ggplot2 | Add geom_text annotations for difference estimates from gcomp results. |
| D-6 | **as_gt() / as_gtsummary() methods for pool objects** | S3 methods that let users pipe directly: `pool_obj |> as_gt()`. Clean API that feels native to the gt/gtsummary ecosystem. | Low | gt, gtsummary | Thin wrappers calling internal tidy -> format -> table pipeline. |

---

## Anti-Features

Features to explicitly NOT build. Common mistakes in this domain that would increase maintenance burden without proportional value.

| # | Anti-Feature | Why Avoid | What to Do Instead |
|---|-------------|-----------|-------------------|
| AF-1 | **Full rtables/tern integration** | rtables has a fundamentally different architecture (cell-based layout engine) from gt/gtsummary (ARD-based pipeline). Supporting both ecosystems doubles the table generation code and testing surface. tern is tightly coupled to rtables. The pharmaverse is converging on ARD+gtsummary for new development (see crane package, gtsummary ARD-first functions). | Commit to the cards/gtsummary/gt stack. Users who need rtables can extract the tidy tibble and build their own rtables layout -- the data is accessible. |
| AF-2 | **Interactive Shiny dashboards** | Already out of scope per PROJECT.md. Teal/teal.modules.clinical covers this space comprehensively. Building interactive features fragments the package's focus and adds massive dependency burden (shiny, teal, plotly). | Export ggplot2 objects that teal modules can consume. Static figures are the deliverable. |
| AF-3 | **Word/PowerPoint/RTF direct export** | gt handles HTML and PDF well. Word/RTF export requires flextable or officer, adding dependencies for a use case better served by separate tools. Regulatory submissions use PDF; internal review uses HTML. | Support gt's existing `gtsave()` for HTML/PDF. Document how to use `gt::as_rtf()` if Word output is needed. Do not add flextable as a dependency. |
| AF-4 | **Safety tables (AE, labs, vital signs)** | Outside rbmi's efficacy analysis domain. Safety reporting is handled by packages like teal.modules.clinical, tern, and gtreg. Building safety tables would blur the package's scope. | Keep scope to efficacy endpoints that flow through rbmi's imputation-analysis-pooling pipeline. |
| AF-5 | **Custom pooling methods** | as_analysis2() hardcodes pooling class assignment based on rbmi method types. Adding custom pooling would require changes to rbmi internals or reimplementation of pool(). | Use rbmi's built-in pooling. If users need custom pooling, they should work with rbmi directly. Document this limitation. |
| AF-6 | **Generic table builder / layout engine** | Building a general-purpose table layout system would compete with gt, rtables, and gtsummary -- all mature, well-funded projects. rbmiUtils should produce tables, not be a table framework. | Depend on gt/gtsummary for table rendering. Provide rbmi-specific content and opinionated defaults, not layout primitives. |
| AF-7 | **Spaghetti/trajectory plots** | Individual patient trajectory plots are a different visualization category from summary-level efficacy plots. They require raw imputed data (not pooled results) and serve a different analytical purpose. | Provide `get_imputed_data()` for users who want to build trajectory plots themselves. Focus figure functions on pooled summary results. |
| AF-8 | **Reimplementing cards/cardx statistical functions** | The temptation exists to build custom ARD creation for ANCOVA, MMRM, etc. But cardx already provides `ard_car_anova()`, `ard_emmeans_contrast()`, and model-based ARD functions. | Use cardx where applicable. Only build custom ARD conversion for rbmi-specific output (pooled MI results, gcomp results) that cardx cannot handle. |

---

## Feature Dependencies

```
                    TS-1: ARD Conversion
                     /          \
                    /            \
            TS-2: Efficacy     D-1: One-call
              Table            efficacy_table()
                |                    |
                v                    v
            TS-6: Column        D-2: Sensitivity
            Formatting          Comparison Table
                                     |
                                     v
                                D-3: MI Metadata
                                   in ARD

            TS-3: Forest Plot
                |
                v
            D-4: Sensitivity Overlay

            TS-5: Responder Bar Chart
                |
                v
            D-5: Annotated Differences

            TS-4: Print/Summary Methods (independent)
```

### Critical Path

1. **TS-4** (print/summary) is independent -- can be built first or in parallel.
2. **TS-1** (ARD) must come before TS-2 and D-1.
3. **TS-3** (forest plot) and **TS-5** (responder bar chart) are independent of the ARD/table path.
4. **TS-6** (column formatting) integrates with TS-2 but builds on existing format_* functions.

### Dependency on Existing Functions

| New Feature | Consumes From Existing Code |
|-------------|----------------------------|
| TS-1 (ARD) | `tidy_pool_obj()` output tibble |
| TS-2 (Efficacy table) | TS-1 ARD output, `extract_trt_effects()`, `extract_lsm()` |
| TS-3 (Forest plot) | `extract_trt_effects()` output |
| TS-4 (Print/summary) | `tidy_pool_obj()`, `format_results()` |
| TS-5 (Responder bar chart) | `gcomp_responder_multi()` pooled results via `tidy_pool_obj()` |
| D-1 (One-call table) | TS-1 + TS-2 + TS-6 combined |
| D-2 (Sensitivity table) | `combine_results()`, TS-1 |

---

## Ecosystem Positioning

### Where rbmiUtils Fits in the Pharmaverse

| Layer | Package(s) | rbmiUtils Role |
|-------|-----------|----------------|
| Imputation | rbmi | Consumer (wraps rbmi::pool output) |
| Analysis Results | cards, cardx | Producer (generates ARD from MI-specific results) |
| Table rendering | gt, gtsummary, crane | Consumer (uses gt/gtsummary to render tables) |
| Figures | ggplot2 | Consumer (uses ggplot2 for forest/bar plots) |
| Interactive | teal, teal.modules.clinical | Out of scope (exports static objects compatible with teal) |
| Safety reporting | tern, gtreg | Out of scope (efficacy only) |

### Competitive Landscape

| Capability | tern + rtables | gtsummary + crane | rbmiUtils (proposed) |
|-----------|---------------|-------------------|---------------------|
| Table engine | rtables (cell-based) | gt (ARD-based) | gt via gtsummary (ARD-based) |
| ARD support | No (predates ARD) | Native | Via cards bridge |
| MI-specific output | Generic (no MI metadata) | Generic | MI-aware (imputations, pooling method, FMI) |
| Efficacy tables | 220+ templates in TLG Catalog | General-purpose | Opinionated for rbmi ANCOVA/MMRM |
| Sensitivity comparison | Manual construction | Manual construction | Built-in (D-2) |
| Forest plots | g_forest() in tern | Not built-in | Built-in for MI visit/subgroup |
| Responder charts | Available in tern | Not built-in | Built-in for gcomp results |
| Setup overhead | High (rtables learning curve) | Medium (ARD pipeline) | Low (one-call from pool object) |

### Key Ecosystem Insight

The pharmaverse is undergoing a transition. The older stack (rtables + tern) provides 220+ table templates but requires learning the rtables cell-based layout system. The newer stack (cards + gtsummary + crane + gt) is ARD-first and gaining adoption rapidly (cards has 52,000+ monthly CRAN downloads). **rbmiUtils should align with the ARD-first stack** because:

1. cards/gtsummary is where Roche, GSK, Novartis, and Pfizer are investing
2. ARD provides traceability and reproducibility mandated by CDISC
3. gtsummary's `tbl_ard_*()` functions accept pre-computed ARDs -- a natural fit for rbmi pooled results
4. crane demonstrates the pattern: extend gtsummary with domain-specific opinionated defaults

---

## Expected Output: Regulatory Efficacy Table

A standard regulatory submission efficacy table (ICH E3 Section 11) for an rbmi analysis looks like this:

```
Table 14.2.1: Primary Efficacy Analysis - Change from Baseline in [Endpoint]
MMRM/ANCOVA with Multiple Imputation (Reference-Based)

                    Placebo        Drug A         Treatment Difference
Visit              LS Mean (SE)   LS Mean (SE)    Est (95% CI)           P-value
---------------------------------------------------------------------------
Week 4             1.20 (0.45)    0.85 (0.43)    -0.35 (-1.57, 0.87)    0.574
Week 8             2.15 (0.52)    1.10 (0.50)    -1.05 (-2.49, 0.39)    0.153
Week 12            3.40 (0.58)    1.75 (0.56)    -1.65 (-3.23, -0.07)   0.041
Week 16            4.10 (0.62)    2.05 (0.60)    -2.05 (-3.74, -0.36)   0.018
Week 24*           5.20 (0.65)    2.75 (0.71)    -2.45 (-4.20, -0.70)   0.006
---------------------------------------------------------------------------
* Primary endpoint

Number of imputations: 100
Pooling method: Rubin's rules
Model: ANCOVA adjusted for baseline, strata, region
Reference-based imputation: Jump-to-reference
```

This table must be producible from a single rbmi pool object. The mapping from existing functions:
- LS Means: `extract_lsm(tidy_result)` per arm
- Treatment Difference: `extract_trt_effects(tidy_result)`
- Formatting: `format_estimate()` + `format_pvalue()`
- Table: gt or gtsummary rendering

---

## MVP Recommendation

For the reporting milestone, prioritize in this order:

### Phase 1: Foundation (must ship)
1. **TS-4: Print/summary methods** -- Independent, low complexity, immediate user value. Improves developer experience for every rbmi user.
2. **TS-1: ARD conversion** -- Enables the entire table pipeline. Without this, no regulatory tables via gtsummary.
3. **TS-3: Forest plot** -- High-visibility deliverable. Independent of ARD path. Users can see it working immediately.

### Phase 2: Tables (core value)
4. **TS-2: Efficacy summary table** -- The flagship deliverable. Consumes TS-1 ARD.
5. **TS-5: Responder bar chart** -- Completes the binary endpoint reporting story alongside gcomp functions.
6. **TS-6: Column formatting** -- Polish for regulatory table output.

### Phase 3: Differentiators (competitive advantage)
7. **D-1: One-call efficacy_table()** -- The "magic" API that makes rbmiUtils feel effortless.
8. **D-2: Sensitivity comparison table** -- Unique to MI workflow, high regulatory value.

### Defer to Post-MVP
- **D-3** (MI metadata in ARD): Nice-to-have, low urgency
- **D-4** (Sensitivity overlay on forest): Advanced visualization, can wait
- **D-5** (Annotated responder chart): Polish feature
- **D-6** (as_gt/as_gtsummary methods): API sugar, can add after core functions work

---

## Technical Decisions Needed Before Implementation

| Decision | Options | Recommendation | Rationale |
|----------|---------|----------------|-----------|
| cards/gtsummary as Imports vs Suggests | Imports (always available) vs Suggests (optional) | **Suggests** for cards, gtsummary, gt | Keeps base rbmiUtils installable without heavy dependencies. Use `rlang::check_installed()` at function entry. |
| ARD function name | `pool_to_ard()` vs `as_ard.pool()` vs `tidy_to_ard()` | **`pool_to_ard()`** | Explicit, not an S3 method (avoids masking cards methods), clearly communicates direction of conversion. |
| Table rendering backend | gt only vs gtsummary | **gtsummary for structure, gt for rendering** | gtsummary provides `tbl_ard_*()` functions purpose-built for ARD consumption. gt handles final rendering. This matches the crane pattern. |
| Forest plot implementation | ggplot2 custom vs forestploter vs ggforestplot | **ggplot2 custom** | ggplot2 is already a Suggests dependency. Custom implementation gives full control over rbmi-specific layout (visits as y-axis, sensitivity overlays). No additional dependency. |
| Responder bar chart | ggplot2 custom vs existing packages | **ggplot2 custom** | Same rationale as forest plot. Simple `geom_col()` + `geom_errorbar()` with treatment arm grouping. |

---

## Complexity Estimates

| Feature | Lines of R Code (est.) | Test Lines (est.) | Documentation | Total Effort |
|---------|----------------------|-------------------|---------------|-------------|
| TS-1: ARD conversion | 80-120 | 60-80 | Vignette section | Medium |
| TS-2: Efficacy table | 150-200 | 80-100 | Vignette + examples | Medium-High |
| TS-3: Forest plot | 100-150 | 40-60 | Examples + vignette section | Medium |
| TS-4: Print/summary | 60-100 | 40-60 | Roxygen only | Low |
| TS-5: Responder bar chart | 80-120 | 40-60 | Examples | Low-Medium |
| TS-6: Column formatting | 40-60 | 30-40 | Roxygen only | Low |
| D-1: One-call table | 60-80 | 40-60 | Vignette highlight | Low-Medium |
| D-2: Sensitivity table | 100-140 | 60-80 | Vignette section | Medium |

---

## Sources

### HIGH Confidence (official docs, codebase inspection)
- rbmiUtils codebase: Direct inspection of R/, tests/, DESCRIPTION
- rbmi pool object structure: [GitHub rbmi/R/pool.R](https://github.com/insightsengineering/rbmi/blob/main/R/pool.R)
- cards package: [GitHub insightsengineering/cards](https://github.com/insightsengineering/cards)
- cards `tidy_as_ard()`: [cards reference docs](https://insightsengineering.github.io/cards/v0.3.0/reference/tidy_as_ard.html)
- gtsummary ARD-first tables: [gtsummary tbl_ard functions](https://www.danieldsjoberg.com/gtsummary/articles/tbl_ard-functions.html)
- cardx package functions: [cardx main docs](https://insightsengineering.github.io/cardx/main/)
- tern TLG catalog: [tern package](https://insightsengineering.github.io/tern/main/)
- crane package: [GitHub insightsengineering/crane](https://github.com/insightsengineering/crane)
- ICH E3 guideline: [FDA E3 Structure and Content of Clinical Study Reports](https://www.fda.gov/media/84857/download)

### MEDIUM Confidence (verified with multiple sources)
- ARD as pharmaverse standard: [R Consortium blog on cards](https://r-consortium.org/posts/supercharging-statistical-analysis-with-ards-and-the-cards-r-package/), [CDISC COSA Spotlight 2025](https://www.danieldsjoberg.com/CDISC-COSA-Spotlight-ARD-gtsummary-2025/slides/)
- cards monthly downloads (52,000+): [R Consortium article](https://r-consortium.org/posts/supercharging-statistical-analysis-with-ards-and-the-cards-r-package/)
- Pharmaverse TLG ecosystem: [pharmaverse.org TLGs](https://pharmaverse.org/e2eclinical/tlg/)
- gtreg for regulatory tables: [gtreg CRAN page](https://cran.r-project.org/web/packages/gtreg/index.html)
- Forest plot packages: [forestploter](https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html), [ggforestplot](https://nightingalehealth.github.io/ggforestplot/index.html)

### LOW Confidence (needs validation during implementation)
- Exact mapping from tidy_pool_obj() columns to ARD structure -- will need prototyping
- Whether gtsummary's tbl_ard_summary() is flexible enough for the regulatory table layout, or whether direct gt construction is needed
- crane's theme pattern applicability to rbmiUtils (crane is Roche-specific; rbmiUtils is open)

---

*Feature landscape research: 2026-02-07*
