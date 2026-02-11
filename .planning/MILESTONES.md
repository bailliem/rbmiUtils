# Project Milestones: rbmiUtils

## v1 Reporting & Robustness (Shipped: 2026-02-08)

**Delivered:** Added a reporting layer bridging rbmi pooled results into publication-ready regulatory tables and forest plots, with hardened foundation functions and pharmaverse ARD integration

**Phases completed:** 1-4 (9 plans total)

**Key accomplishments:**
- Hardened parameter parsing, class detection, gcomp validation, and storage round-trip integrity across the foundation layer
- Added informative print/summary S3 methods for pool and analysis objects with cli formatting and significance flags
- Built pool_to_ard() converting rbmi results to pharmaverse ARD standard via cards package
- Created efficacy_table() producing regulatory-style gt tables with LS means, treatment differences, CIs, and p-values by visit
- Implemented plot_forest() with three-panel patchwork composition (table | forest | p-values), treatment and LSM modes, and Okabe-Ito colorblind palette

**Stats:**
- 67 files created/modified
- 3,976 lines of R source code, 4,916 lines of test code
- 4 phases, 9 plans, 48 commits
- 2 days (2026-02-07 to 2026-02-08)

**Git range:** `9088fbd` → `d197e29`

**What's next:** v2 features — responder bar charts, sensitivity analysis overlays, column formatting controls, MI-specific ARD metadata, draws/imputation print helpers

---

## v2 Documentation & Hardening (Shipped: 2026-02-10)

**Delivered:** Made rbmiUtils discoverable and trustworthy with a polished pkgdown site, end-to-end vignette, hardened data prep functions, and comprehensive documentation

**Phases completed:** 5-7 (7 plans total)

**Key accomplishments:**
- Hardened validate_data() and prepare_data_ice() with cli messaging and 6 new validation checks (malformed interaction terms, NULL strategy, character visit warnings, empty data frames, all-NA covariates, batched coercion warnings)
- Added 13 edge case tests covering single subject, single visit, all-NA outcome, and all-complete data scenarios (95 total data prep tests)
- Created end-to-end pipeline vignette from ADEFF data through rbmi draws/impute/analyse/pool to efficacy_table() and plot_forest(), with binary responder appendix using beeca
- Pre-rendered documentation images for README visual teasers and roxygen help page examples
- Restructured NEWS.md into versioned changelog (0.2.0/0.1.0) with inline rbmi/beeca cross-references across all vignettes
- Configured complete pkgdown site with hex logo, 9-group reference index, navbar, open graph social cards, and pharmaverse footer

**Stats:**
- 69 files created/modified
- 4,146 lines of R source code, 5,285 lines of test code
- 3 phases, 7 plans
- 2 days (2026-02-08 to 2026-02-10)

**Git range:** `01c2d22` → `51c4348`

**What's next:** v3 features — responder bar charts, sensitivity analysis overlays, column formatting controls, MI-specific ARD metadata, draws/imputation helpers

---


## v3 ARD Enrichment & Polish (Shipped: 2026-02-11)

**Delivered:** Enriched ARD output with MI diagnostic metadata, added pipeline introspection helpers, polished tables and plots to publication quality, and overhauled all documentation with realistic examples

**Phases completed:** 8-11 (9 plans total)

**Key accomplishments:**
- Enriched pool_to_ard() with 7 MI diagnostic statistics per parameter (FMI, lambda, RIV, Barnard-Rubin df, relative efficiency) via compute_rubin_diagnostics(), with curated ARD rows passing cards validation
- Added describe_draws() and describe_imputation() introspection helpers with cli-formatted print methods, MCMC convergence diagnostics (ESS, Rhat), and missingness breakdown by visit/arm
- Added publication styling controls to efficacy_table() (font_family, font_size, row_padding) and plot_forest() (font_family, panel_widths) with backward-compatible NULL defaults
- Overhauled README with complete ADEFF-through-pipeline Quick Start, upgraded 5 function examples to use realistic package datasets (ADMI/ADEFF)
- Created MI diagnostics and pipeline inspection vignette, regenerated forest plot images, bumped version to v0.3.0 with comprehensive NEWS.md entries
- Added 400+ new test expectations across 9 plans (98 data prep tests total, 6,668 total test lines) with zero R CMD check errors/warnings

**Stats:**
- 57 files changed, 9,700 insertions, 145 deletions
- 4,980 lines of R source code, 6,668 lines of test code
- 4 phases, 9 plans, 48 commits
- 2 days (2026-02-10 to 2026-02-11)

**Git range:** `59e4ff4` → `5d8b503`

**What's next:** v4 features — responder bar charts, sensitivity analysis overlays, column formatting controls, as_gt()/as_gtsummary() S3 methods, BMLMI diagnostics

---

