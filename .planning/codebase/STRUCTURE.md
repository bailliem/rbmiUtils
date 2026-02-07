# Codebase Structure

**Analysis Date:** 2026-02-07

## Directory Layout

```
rbmiUtils-gsd/
├── R/                    # Core package code (11 files, ~2612 lines)
├── tests/                # Test suite
│   ├── testthat/        # Unit and integration tests
│   └── testthat.R       # Test configuration
├── vignettes/           # Package documentation and examples
├── data/                # Example datasets (ADEFF, ADMI)
├── data-raw/            # Scripts to generate example datasets
├── inst/                # Installation files
├── man/                 # Auto-generated documentation (Roxygen)
├── docs/                # Pkgdown website
├── DESCRIPTION          # Package metadata
├── NAMESPACE            # Package exports
├── README.md            # User-facing documentation
└── .planning/codebase/  # GSD codebase mapping documents
```

## Directory Purposes

**R/ - Core Implementation:**
- Purpose: All source code for user-facing functions and helpers
- Contains: 11 R files implementing data validation, analysis, result processing, and utilities
- Key files: `analyse_mi_data.R` (442 lines), `data_helpers.R` (553 lines), `result_helpers.R` (386 lines)

**tests/testthat/ - Test Suite:**
- Purpose: Unit and integration tests for all public functions
- Contains: 9 test files covering data validation, analysis functions, results, formatting, storage
- Test files: Named `test-{module}.R` matching source file names (e.g., `test-analyse_mi_data.R`)
- Examples: `test-analyse_mi_data.R`, `test-data_helpers.R`, `test-formatting.R`, `test-imputation_storage.R`

**vignettes/ - User Guides:**
- Purpose: High-level workflow documentation with examples
- Contains: 3 R Markdown vignettes
  - `data-preparation.Rmd` - data validation and prep workflow
  - `analyse2.Rmd` - analysis execution and pooling workflow
  - `efficient-storage.Rmd` - reduce/expand storage optimization

**data/ - Example Datasets:**
- Purpose: Bundled datasets for examples and vignettes
- Contains: ADEFF (original efficacy data with missingness), ADMI (100,000-row fully imputed dataset)
- Use: Load with `data("ADMI")` in examples and vignettes

**man/ - Documentation:**
- Purpose: Auto-generated Roxygen documentation
- Contains: .Rd files corresponding to exported functions (not edited directly)
- Generated from: Roxygen comments in R/ files

## Key File Locations

**Entry Points:**
- `R/analyse_mi_data.R` (line 79): `analyse_mi_data()` - main analysis workflow function
- `R/data_helpers.R` (line 64): `validate_data()` - pre-flight data validation
- `R/utils.R` (line 52): `get_imputed_data()` - extract imputed datasets from rbmi objects
- `R/tidiers.R` (line 88): `tidy_pool_obj()` - format pooled results for publication

**Configuration:**
- `DESCRIPTION` - package metadata, dependencies (rbmi, dplyr, assertthat, beeca, rlang, tidyr)
- `NAMESPACE` - exported functions (23 exports defined)
- `_pkgdown.yml` - website configuration

**Core Logic by Module:**

**Data Preparation/Validation:**
- `R/data_helpers.R` (553 lines)
  - `validate_data()` - comprehensive pre-imputation validation
  - `prepare_data_ice()` - intercurrent event data construction
  - `summarise_missingness()` - missing data pattern reporting

**Analysis Execution:**
- `R/analyse_mi_data.R` (442 lines)
  - `analyse_mi_data()` - apply analysis function to each imputation
  - `as_analysis2()` - construct analysis object
  - `print.analysis()`, `summary.analysis()` - S3 methods

**Analysis Functions (Binary Outcomes):**
- `R/analysis_utils.R` (229 lines)
  - `gcomp_responder()` - g-computation for single visit binary analysis
  - `gcomp_responder_multi()` - multi-visit wrapper
  - Helper functions: `as_simple_formula2()`, `extract_covariates2()`

**Result Processing:**
- `R/tidiers.R` (167 lines): `tidy_pool_obj()` - transform pooled object to tidy tibble
- `R/formatting.R` (261 lines):
  - `format_pvalue()` - p-value formatting (e.g., "< 0.001")
  - `format_estimate()` - estimate formatting with specified decimals
  - `format_results_table()` - full table formatting
- `R/result_helpers.R` (386 lines):
  - `create_impid()` - add IMPID to list of datasets
  - `combine_results()` - combine multiple endpoint results
  - `format_results()` - apply formatting to tidy results

**Storage Optimization:**
- `R/imputation_storage.R` (327 lines)
  - `reduce_imputed_data()` - compress by keeping only originally missing rows
  - `expand_imputed_data()` - reconstruct full imputed dataset from reduced form

**Utilities:**
- `R/utils.R` (174 lines): `get_imputed_data()`, `gcomp_binary()`
- `R/formatting.R`: Formatting helpers (format_pvalue, format_estimate, format_results_table)

**Package Metadata:**
- `R/rbmiUtils-package.R` (31 lines): globalVariables declaration, package documentation
- `R/ADMI.R` (22 lines): ADMI dataset documentation
- `R/ADEFF.R` (20 lines): ADEFF dataset documentation

## Naming Conventions

**Files:**
- Source code: `{function_group}.R` (e.g., `data_helpers.R`, `analyse_mi_data.R`)
- Tests: `test-{function_group}.R` matching source modules
- Vignettes: `{workflow-name}.Rmd` (e.g., `data-preparation.Rmd`, `analyse2.Rmd`)
- Example datasets: `{DATASET}.R` documenting structure (e.g., `ADMI.R`, `ADEFF.R`)

**Functions:**
- Snake_case for multi-word functions: `analyse_mi_data()`, `prepare_data_ice()`, `reduce_imputed_data()`
- Descriptive action verbs: `validate_`, `prepare_`, `get_`, `format_`, `tidy_`, `expand_`, `reduce_`, `combine_`, `summarise_`
- S3 methods: `print.analysis()`, `summary.analysis()`

**Variables (within functions):**
- `data` - input dataset
- `vars` - vars object from rbmi::set_vars()
- `method` - imputation method object from rbmi
- `results` - analysis results (list or tibble)
- `IMPID` - imputation identifier column (uppercase, special meaning)
- `outcome`, `visit`, `group`, `subjid` - variable specifications from vars object
- Standard abbreviations: `df` (data frame), `col` (column), `msg` (message)

**Constants/Special Columns:**
- `IMPID` - imputation ID column (must be present for analyse_mi_data)
- `original_id` - temporary column mapping internal to external subject IDs
- `internal_id` - temporary column for internal subject IDs in get_imputed_data
- `delta` - sensitivity analysis adjustment column (optional)

## Where to Add New Code

**New Analysis Function (similar to gcomp_responder):**
- Implementation: `R/analysis_utils.R` (if binary/responder), or new `R/{function_group}.R`
- Tests: `tests/testthat/test-analysis_utils.R` or `tests/testthat/test-{function_group}.R`
- Pattern: Function takes (data, vars, ...), returns named list with $trt, optionally other effects
- Integration: Pass to `analyse_mi_data(..., fun = your_function)`
- Documentation: Roxygen @export, @param, @return, @examples

**New Data Validation Check:**
- Location: `R/data_helpers.R` within `validate_data()` function
- Pattern: Append to `issues` character vector, let function report all at once
- Test: `tests/testthat/test-data_helpers.R`

**New Formatting Function:**
- Location: `R/formatting.R` (publication formatting) or `R/result_helpers.R` (result combination)
- Pattern: Take tibble, return tibble with modified columns
- Export: Add `@export` to Roxygen
- Test: `tests/testthat/test-formatting.R` or `tests/testthat/test-result_helpers.R`

**Utilities/Helpers:**
- Location: `R/utils.R` or appropriate module file
- Pattern: Standalone functions or internal helpers (no @export)
- Test: Place in test file matching module

**Example/Documentation:**
- Vignettes: Add `{workflow}.Rmd` to `vignettes/`
- Build site: `pkgdown::build_site()`
- Include in examples: Use `@examples` in Roxygen documentation

## Special Directories

**man/ - Generated Documentation:**
- Purpose: Contains .Rd files auto-generated from Roxygen comments
- Generated: Yes, from source code via roxygen2
- Committed: Yes, to support offline documentation
- Do not edit directly: Edit Roxygen comments in R/ files instead

**docs/ - Website:**
- Purpose: Pkgdown-generated website
- Generated: Yes, via pkgdown::build_site()
- Committed: Yes, for GitHub Pages hosting
- Contents: Reference docs, articles (vignettes), news

**.planning/codebase/ - GSD Codebase Maps:**
- Purpose: Generated by GSD system to document architecture, structure, conventions, testing patterns
- Generated: Yes, by /gsd:map-codebase command
- Committed: Yes, for context in future GSD phases

**data/ - Data Objects:**
- Purpose: Compiled binary R datasets (.rda)
- Generated: Yes, from data-raw/ scripts
- Committed: Yes, to support data("ADMI") loading
- Contents: ADMI, ADEFF example datasets

**data-raw/ - Data Generation Scripts:**
- Purpose: Scripts to generate datasets in data/
- Contains: R scripts that create ADMI, ADEFF
- Run with: devtools::load_all() or R CMD build

---

*Structure analysis: 2026-02-07*
