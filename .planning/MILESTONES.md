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
