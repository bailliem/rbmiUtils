# Project Research Summary

**Project:** rbmiUtils reporting milestone (ARD/gtsummary integration)
**Domain:** Clinical trial reporting utilities for reference-based multiple imputation analysis
**Researched:** 2026-02-07
**Confidence:** MEDIUM-HIGH

## Executive Summary

rbmiUtils bridges rbmi's multiple imputation analysis pipeline to regulatory-quality reporting. The research confirms that the package should adopt the pharmaverse's ARD-first workflow using cards/gtsummary/gt for tables and ggplot2 for figures. This aligns with where Roche, GSK, Novartis, and Pfizer are investing and provides the traceability/reproducibility mandated by CDISC standards.

The recommended architecture adds three new layers (ARD Conversion, Table Generation, Visualization) on top of six existing layers, with a clean boundary at the tidy tibble produced by `tidy_pool_obj()`. The core technical challenge is converting rbmi's pooled results (wide format: est, se, lci, uci, pval) into the cards ARD format (long format: stat_name/stat/stat_label with list-columns). The conversion uses `cards::ard_identity()` for pre-calculated statistics, pivoting each row of tidy pool output into 5 ARD rows (one per statistic).

Critical risks center on fragile parameter parsing with underscored visit names, rbmi class hierarchy coupling via `class(method)[[2]]`, and ARD column contract violations. The recommended mitigation is to harden the existing codebase first (fix parsing with regex, wrap `rbmi::analyse()` to eliminate class indexing), then build ARD conversion and reporting layers on stable foundations. Phase 1 must address these structural issues before Phase 2 adds cards/gtsummary dependencies.

## Key Findings

### Recommended Stack

The research identifies cards (v0.7.1) as the foundation for ARD conversion, with gtsummary (v2.5.0) and gt (v1.3.0) providing table generation and multi-format rendering. All three packages are actively maintained with pharma-specific features and official CDISC backing. cards has 52,000+ monthly CRAN downloads and co-development by Roche/GSK/Novartis.

**Core technologies:**
- **cards (>= 0.7.0)** — ARD data structure and `ard_identity()` for wrapping pre-calculated statistics. Only hard dependency needed; minimal transitive dependencies (dplyr, tidyr, rlang, cli).
- **gtsummary (>= 2.0.0)** — ARD-first table functions (`tbl_ard_wide_summary()` for efficacy tables). Suggests dependency to keep base install light.
- **gt (>= 1.0.0)** — Multi-format rendering (HTML, RTF, LaTeX/PDF, Word). Stable APIs for regulatory submissions. Also Suggests.
- **ggplot2 (>= 3.5.0)** — Already in Suggests. Sufficient for forest plots and responder bar charts without specialized packages.
- **forestploter (>= 1.0.0)** — Optional Suggests for advanced forest plots with table-aligned layout. CRAN maintained (April 2025).

**Avoid:**
- rtables/tern (different ecosystem, incompatible with ARD paradigm)
- tfrmt/chevron (adds complexity without clear benefit when gtsummary already consumes ARDs)
- flextable/r2rtf as hard dependencies (gt handles primary use cases; fall back only when needed)

### Expected Features

**Must have (table stakes):**
- **TS-1: ARD conversion from tidy pool results** — The pharmaverse has standardized on ARD as the interchange format. Any package producing clinical trial results needs ARD output to be interoperable. Maps naturally from tidy tibble (est, se, lci, uci, pval) to ARD (stat_name/stat structure).
- **TS-2: Efficacy summary table (regulatory style)** — Every clinical study report includes a primary efficacy table showing LS means by arm, treatment differences, CIs, and p-values by visit. Standard layout: rows = visits, columns = Reference LS Mean, Treatment LS Mean, Difference (95% CI), P-value.
- **TS-3: Forest plot for treatment effects** — Standard visualization showing treatment effects across visits with point estimates and CIs. Expected in every regulatory submission with longitudinal efficacy data.
- **TS-4: Print/summary methods for pool objects** — Users expect enhanced print methods showing rounded estimates, formatted CIs, and clear parameter labeling at the console.
- **TS-5: Responder analysis bar chart** — Binary endpoint reporting requires grouped bar charts showing proportion responding by treatment arm across visits. Visual counterpart to `gcomp_responder_multi()` results.
- **TS-6: Column formatting controls for tables** — Regulatory tables require specific decimal precision, CI bracket styles, p-value thresholds. Extend existing formatters to integrate with gt/gtsummary theming.

**Should have (competitive advantage):**
- **D-1: One-call efficacy table from pool object** — `efficacy_table(pool_obj)` that goes directly from `rbmi::pool()` to gt table with no intermediate steps. Competing workflows require 5-10 lines. A function that "just works" for rbmi is a strong differentiator.
- **D-2: Sensitivity analysis comparison table** — Side-by-side table comparing primary vs. tipping point / delta-adjusted results. Unique to MI workflow where sensitivity analyses are mandatory. No existing package does this for rbmi output.
- **D-3: Integrated ARD with MI metadata** — Extend ARD with MI-specific metadata (number of imputations, pooling method, fraction of missing information) for self-documenting regulatory review.

**Defer (v2+):**
- Interactive Shiny dashboards (teal/teal.modules.clinical covers this)
- Safety tables (AE, labs, vital signs — outside rbmi efficacy domain)
- Word/PowerPoint direct export beyond gt's built-in support
- Spaghetti/trajectory plots (different visualization category from pooled results)

### Architecture Approach

The clean architecture boundary sits at the tidy tibble produced by `tidy_pool_obj()`. Analysis functions (layers 1-4) should never know about cards/gtsummary. Reporting functions (layers 5-7) should never reach back into rbmi imputation mechanics. This enables testing reporting functions with mocked tidy tibbles without needing the full rbmi pipeline.

**Major components:**
1. **ARD Conversion Layer (new)** — Converts tidy tibbles to ARD format using `cards::ard_identity()`. Each row of tidy pool output (one parameter with 5 statistics) becomes 5 rows of ARD (one per statistic with stat_name, stat_label, stat as list-column). Depends on cards (Suggests).
2. **Table Generation Layer (new)** — Produces regulatory tables from ARDs using `gtsummary::tbl_ard_wide_summary()` for efficacy tables (treatment effects by visit) and `tbl_ard_summary()` for responder tables. Renders via gt to HTML/RTF/PDF/Word. Depends on gtsummary + gt (Suggests).
3. **Visualization Layer (new)** — Produces forest plots (`plot_forest()`) and responder bar charts (`plot_responder_bar()`) using ggplot2. Works directly from tidy tibbles, independent of ARD path. Depends on ggplot2 (already Suggests).

**Key pattern:** Tidy-first, ARD-second. Always convert pool objects to tidy tibbles first, then to ARD. Never bypass the tidy tibble contract. Use soft dependency checks (`rlang::check_installed()`) at function entry for all Suggests packages.

### Critical Pitfalls

1. **ARD Column Contract Violations** — The `stat` column must be a list-column, not numeric vector. Missing `context`, `fmt_fn`, or wrong column types cause `tbl_ard_summary()` to fail silently or error cryptically. Prevention: Use `cards::as_card()` to guarantee structure, validate output class, map statistics explicitly (est -> "estimate", se -> "std.error", lci -> "conf.low", uci -> "conf.high", pval -> "p.value").

2. **Parameter Name Parsing Breaks on Underscored Visits** — `tidy_pool_obj()` uses `tidyr::separate(parameter, sep = "_")` which fails when visit names contain underscores (e.g., "Week_24", "Follow_Up_Visit"). The delimiter collision causes silent data corruption where `parameter_type` gets partial strings and `visit` is truncated. Prevention: Replace with structured regex parsing (`"^(trt|lsm)_(ref_|alt_)?(.+)$"`) or migrate to `tidyr::separate_wider_regex()` with named capture groups.

3. **rbmi Class Hierarchy Indexing with `class(method)[[2]]`** — `analyse_mi_data()` and `as_analysis2()` use positional indexing to determine method type (bayes, approxbayes, condmean, bmlmi). If rbmi changes class hierarchy order, this silently assigns wrong pooling methods, corrupting treatment effects and p-values. Prevention: Replace with `inherits()` checks or wrap `rbmi::analyse()` directly to avoid reimplementing the mapping.

4. **beeca Output Format Coupling** — `gcomp_responder()` directly accesses columns `STAT`, `STATVAL`, `TRTVAL` from `beeca::get_marginal_effect()$marginal_results` with no version constraint. Column renames break the entire analysis pipeline. Prevention: Pin beeca version, validate output schema immediately after beeca call, abstract column mapping to single location.

5. **gtsummary Template Brittleness** — Regulatory efficacy tables have non-standard layouts (multiple statistics per cell, visit grouping, specific headers/footers). gtsummary excels at Table 1 but complex layouts may require heavy `modify_*()` customization or direct gt construction. Prevention: Prototype exact target layout before coding, use `tbl_ard_wide_summary()` for multi-statistic columns, consider direct gt assembly for complex regulatory formats.

## Implications for Roadmap

Based on research, the roadmap should prioritize hardening fragile foundations before adding new reporting dependencies. Critical structural issues (parameter parsing, class indexing, beeca coupling) will propagate into ARD conversion if not fixed first.

### Phase 1: Foundation Hardening
**Rationale:** Fix fragile parameter parsing, class hierarchy indexing, and beeca coupling before building on top. These issues are observable in existing code and will corrupt ARD conversion if not addressed.
**Delivers:** Stable tidy tibble contract, robust parameter parsing with regex, wrapped `rbmi::analyse()` eliminating class indexing, validated beeca output, delta uniqueness checks.
**Addresses:** No new features — this is technical debt paydown.
**Avoids:** Pitfall 2 (parameter parsing), Pitfall 3 (class indexing), Pitfall 4 (beeca coupling), Pitfall 10 (delta validation), Pitfall 11 (string key collisions).
**Research flag:** Standard refactoring patterns. No deep research needed. Leverage existing test suite.

### Phase 2: Enhanced Print/Summary Methods
**Rationale:** Independent of reporting stack. Improves developer experience immediately. Low risk, high value. Can run in parallel with or before Phase 1.
**Delivers:** Enhanced `print.analysis()` and `summary.analysis()` showing parameter preview, custom `"tidy_pool"` class with print method, descriptive helpers (`describe_pool()`, `describe_draws()`, `describe_imputation()`).
**Addresses:** TS-4 (print/summary methods).
**Avoids:** Pitfall 7 (S3 dispatch conflicts — use wrapper class for rbmi-owned classes).
**Research flag:** No research needed. Standard S3 method patterns.

### Phase 3: ARD Conversion Layer
**Rationale:** Foundation for all table generation. Must come before tables. Depends on Phase 1 providing stable tidy tibbles.
**Delivers:** `pool_to_ard()` and `tidy_to_ard()` functions converting rbmi results to cards ARD format. Validates ARD structure, maps statistics explicitly, includes MI metadata.
**Addresses:** TS-1 (ARD conversion), D-3 (MI metadata in ARD).
**Uses:** cards (Suggests).
**Implements:** ARD Conversion Layer from ARCHITECTURE.md.
**Avoids:** Pitfall 1 (ARD column contract violations — use `as_card()`, validate output class).
**Research flag:** May need deeper research into `cards::as_card()` vs `cards::tidy_as_ard()` vs `cards::ard_identity()` selection. API stability across cards v0.5-v0.7 should be verified.

### Phase 4: Table Generation
**Rationale:** Flagship deliverable. Consumes ARDs from Phase 3. Can run in parallel with Phase 5 (figures).
**Delivers:** `tbl_efficacy()` for regulatory efficacy summary tables, `tbl_responder()` for responder analysis tables. Column formatting integration (TS-6). One-call convenience function `efficacy_table()` (D-1).
**Addresses:** TS-2 (efficacy table), TS-6 (column formatting), D-1 (one-call table), D-2 (sensitivity comparison table).
**Uses:** gtsummary, gt (Suggests).
**Implements:** Table Generation Layer from ARCHITECTURE.md.
**Avoids:** Pitfall 5 (gtsummary template brittleness — prototype layout first, consider direct gt for complex formats), Pitfall 8 (Suggests breaks CRAN — guard all calls, run `R CMD check --no-suggests` in CI).
**Research flag:** Needs deeper research during phase planning. Which `tbl_ard_*()` function fits efficacy tables? Does `tbl_ard_wide_summary()` handle mixed continuous/treatment-effect rows or is custom gt assembly needed? Flag for `/gsd:research-phase`.

### Phase 5: Visualization
**Rationale:** Can run in parallel with Phase 4. Uses tidy tibbles directly (not ARDs). Independent of table generation.
**Delivers:** `plot_forest()` for treatment effect forest plots, `plot_responder_bar()` for responder analysis bar charts. Return ggplot objects for user customization.
**Addresses:** TS-3 (forest plot), TS-5 (responder bar chart).
**Uses:** ggplot2 (already Suggests).
**Implements:** Visualization Layer from ARCHITECTURE.md.
**Avoids:** Pitfall 6 (forest plot scale/ordering/reference line — parameterize scale, always reverse y-axis, return ggplot object), Pitfall 12 (responder bar chart small denominators — annotate n/N, add CIs).
**Research flag:** Standard ggplot2 patterns. No deep research needed.

### Phase Ordering Rationale

- **Phase 1 first because** structural fragility will propagate. Parameter parsing (Pitfall 2) and class indexing (Pitfall 3) corruption will appear in ARD conversion if not fixed first. Wrapping `rbmi::analyse()` eliminates an entire category of maintenance burden.
- **Phase 2 independent** — can run anytime, provides immediate value. Low risk of blocking other work.
- **Phase 3 before Phase 4** — Tables consume ARDs. ARD conversion must work before table generation can start.
- **Phase 4 and Phase 5 parallel** — Figures use tidy tibbles directly, not ARDs. No dependency between tables and figures.
- **Defer advanced features (D-4, D-5, D-6)** to post-MVP. Sensitivity overlay on forest plots, annotated responder charts, and `as_gt()` S3 methods are polish features that can wait.

### Research Flags

Phases likely needing deeper research during planning:
- **Phase 3 (ARD Conversion):** Cards API selection (`as_card()` vs `tidy_as_ard()` vs `ard_identity()`) and stability across versions. Exact mapping from tidy pool columns to ARD structure needs prototyping.
- **Phase 4 (Table Generation):** gtsummary function selection and template layout. Whether `tbl_ard_wide_summary()` is flexible enough for regulatory tables or whether direct gt construction is needed. Likely needs `/gsd:research-phase` during planning.

Phases with standard patterns (skip research-phase):
- **Phase 1 (Hardening):** Standard refactoring. Regex parsing, S3 method best practices, input validation are well-documented patterns.
- **Phase 2 (Print/Summary):** Standard S3 method patterns. No novel techniques.
- **Phase 5 (Visualization):** Standard ggplot2 patterns for forest plots and bar charts. Clinical trial visualization conventions are well-established.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | CRAN verified versions (cards 0.7.1, gtsummary 2.5.0, gt 1.3.0). Official docs confirm ARD workflow. Pharmaverse backing from Roche/GSK/Novartis. |
| Features | MEDIUM-HIGH | Table stakes features (TS-1 through TS-6) validated against codebase needs and ecosystem standards. Differentiators (D-1, D-2, D-3) are rbmi-specific and have no direct precedent but are logical extensions. |
| Architecture | HIGH | Clean boundary at tidy tibble verified from existing code. ARD column structure confirmed from cards official docs. gtsummary ARD consumption patterns verified from vignettes and workshop materials. |
| Pitfalls | MEDIUM-HIGH | Critical pitfalls (1-4) directly observable in codebase with known fragility patterns. Moderate pitfalls (5-9) grounded in ecosystem documentation and domain knowledge. Minor pitfalls (10-13) confirmed via code inspection. |

**Overall confidence:** MEDIUM-HIGH

The stack recommendations (cards/gtsummary/gt) are HIGH confidence — these are official, CRAN-stable packages with clear documentation. The architecture approach (tidy tibble boundary, ARD conversion layer) is HIGH confidence — the data flow is well-defined and matches ecosystem patterns. The feature prioritization is MEDIUM-HIGH — table stakes features are validated against ecosystem expectations, but exact gtsummary template implementation will need iteration. The pitfalls are MEDIUM-HIGH — structural fragilities are directly observable in code, ecosystem-level pitfalls are grounded in documentation, but severity estimates for some (e.g., gtsummary template brittleness) are inferred rather than empirically tested.

### Gaps to Address

- **Exact ARD conversion implementation:** Which cards function (`as_card()`, `tidy_as_ard()`, `ard_identity()`) is best for rbmi pooled results? Needs prototyping during Phase 3 planning.
- **gtsummary template flexibility:** Does `tbl_ard_wide_summary()` handle the specific regulatory table layout (LS means + treatment effects in same table, multiple statistics per cell, visit grouping)? May need fallback to direct gt construction. Validate during Phase 4 planning with `/gsd:research-phase`.
- **forestploter vs ggplot2 for forest plots:** Research identified forestploter (1.1.3) for table-aligned layouts, but ggplot2 custom implementation gives full control and avoids dependency. Decision deferred to Phase 5 planning.
- **rbmi::as_analysis() export status:** Is `as_analysis()` exported publicly or internal-only? Determines whether Phase 1 can use it directly or must maintain `as_analysis2()` copy. Check during Phase 1 start.
- **cards API stability:** Is `as_card()` stable across cards v0.5 to v0.7? May need version guards. Validate during Phase 3 start.

## Sources

### Primary (HIGH confidence)
- [cards 0.7.1 on CRAN](https://cran.r-project.org/web/packages/cards/index.html) — Version verification, ARD structure
- [cards package docs](https://insightsengineering.github.io/cards/main/index.html) — `ard_identity()`, `as_card()`, ARD column contract
- [gtsummary 2.5.0 on CRAN](https://cran.r-project.org/web/packages/gtsummary/index.html) — Version verification
- [gtsummary ARD-first tables](https://www.danieldsjoberg.com/gtsummary/articles/tbl_ard-functions.html) — `tbl_ard_wide_summary()`, ARD consumption
- [gt 1.3.0 on CRAN](https://cran.r-project.org/web/packages/gt/index.html) — Multi-format rendering
- [rbmi quickstart vignette](https://cran.r-project.org/web/packages/rbmi/vignettes/quickstart.html) — `analyse()` and `pool()` workflow
- [rbmi pool.R source](https://github.com/insightsengineering/rbmi/blob/main/R/pool.R) — Pool object structure
- Codebase analysis: `/Users/bailliem/R-projects/rbmiUtils-gsd/R/tidiers.R`, `R/analyse_mi_data.R`, `R/analysis_utils.R`, `R/formatting.R`, `.planning/codebase/CONCERNS.md`

### Secondary (MEDIUM confidence)
- [CDISC COSA Spotlight: ARD + gtsummary 2025](https://www.danieldsjoberg.com/CDISC-COSA-Spotlight-ARD-gtsummary-2025/) — End-to-end ARD workflow, ecosystem adoption
- [ARD PHUSE Workshop 2025](https://www.danieldsjoberg.com/ARD-PHUSE-workshop-2025/) — cards overview, cardx patterns
- [R Consortium: Supercharging Statistical Analysis with ARDs](https://r-consortium.org/posts/supercharging-statistical-analysis-with-ards-and-the-cards-r-package/) — Ecosystem overview, adoption data (52k+ monthly downloads)
- [pharmaverse.org TLG packages](https://pharmaverse.org/e2eclinical/tlg/) — Ecosystem positioning
- [forestploter 1.1.3 on CRAN](https://cran.r-project.org/web/packages/forestploter/index.html) — Forest plot implementation
- [Advanced R: S3](https://adv-r.hadley.nz/s3.html) — S3 dispatch, method conflicts
- [R Packages (2e) Dependencies](https://r-pkgs.org/dependencies-in-practice.html) — Imports vs Suggests, CRAN policy

### Tertiary (LOW confidence, needs validation)
- Exact mapping from tidy_pool_obj() columns to ARD structure — will need prototyping
- Whether gtsummary's `tbl_ard_wide_summary()` handles specific regulatory table layout — needs testing during Phase 4
- cards API stability across v0.5-v0.7 — inferred stable but not directly tested

---
*Research completed: 2026-02-07*
*Ready for roadmap: yes*
