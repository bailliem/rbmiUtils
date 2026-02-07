# Architecture

**Analysis Date:** 2026-02-07

## Pattern Overview

**Overall:** Functional, modular utility library extending rbmi for clinical trial workflows.

**Key Characteristics:**
- Layered architecture separating data preparation, analysis execution, results pooling, and output formatting
- Wrapper functions around rbmi core functions to simplify common workflows
- Heavy use of data validation and assertion patterns (assertthat package)
- Functional composition with higher-order functions (apply analysis function to each imputation)

## Layers

**Data Preparation Layer:**
- Purpose: Validate and prepare data before imputation
- Location: `R/data_helpers.R`
- Contains: Data validation, ICE (intercurrent event) preparation, data summarization
- Depends on: rbmi::set_vars(), dplyr, tidyr
- Used by: Analysis workflow entry points

**Analysis Execution Layer:**
- Purpose: Apply user-defined analysis functions across multiple imputations
- Location: `R/analyse_mi_data.R`
- Contains: Main `analyse_mi_data()` function, analysis object construction, S3 methods for print/summary
- Depends on: Data preparation layer, user-supplied analysis functions (e.g., rbmi::ancova)
- Used by: Results pooling, result helpers

**Analysis Utilities Layer:**
- Purpose: Specialized analysis functions for specific use cases (binary outcomes, g-computation)
- Location: `R/analysis_utils.R`
- Contains: `gcomp_responder()`, `gcomp_responder_multi()`, `gcomp_binary()` - binary outcome analysis
- Depends on: beeca package for marginal effect estimation
- Used by: As alternative analysis functions passed to `analyse_mi_data()`

**Results Processing Layer:**
- Purpose: Transform rbmi pooled objects into tidy, publication-ready formats
- Location: `R/tidiers.R`, `R/formatting.R`, `R/result_helpers.R`
- Contains: Tidying, formatting, result combination and annotation
- Depends on: dplyr, tidyr, rbmi pooled objects
- Used by: Final reporting stage

**Storage Management Layer:**
- Purpose: Optimize storage of imputed datasets through reduction and expansion
- Location: `R/imputation_storage.R`
- Contains: `reduce_imputed_data()` and `expand_imputed_data()` for efficient MI storage
- Depends on: Data preparation layer, rbmi imputation objects
- Used by: Workflows requiring persistent storage of large MI datasets

**Utility Functions Layer:**
- Purpose: Helper functions for common tasks
- Location: `R/utils.R`, `R/result_helpers.R`
- Contains: Data extraction, ID mapping, IMPID creation
- Depends on: rbmi, dplyr, purrr
- Used by: Data preparation and analysis layers

## Data Flow

**Standard Analysis Workflow:**

1. **Input:** Raw clinical trial dataset (ADEFF-like structure with missing values)
2. **Validation:** `validate_data()` checks structure, types, missing patterns
3. **ICE Preparation (optional):** `prepare_data_ice()` creates intercurrent event data
4. **Imputation:** User calls `rbmi::draws()` and `rbmi::impute()` to generate multiple imputations
5. **Data Extraction:** `get_imputed_data()` extracts imputed datasets with original subject IDs mapped
6. **Analysis:** `analyse_mi_data()` applies analysis function (e.g., ancova) to each imputation
7. **Pooling:** `rbmi::pool()` pools results across imputations (Rubin's rules or other methods)
8. **Tidying:** `tidy_pool_obj()` transforms pooled object to publication-ready tibble
9. **Formatting:** `format_results()` applies final formatting (decimals, CI format, p-values)
10. **Combination (optional):** `combine_results()` combines multiple endpoint/subgroup analyses

**Storage Optimization Workflow (alternative):**

1. Extract imputed data from rbmi imputation object
2. `reduce_imputed_data()` - compress by keeping only originally missing rows per imputation
3. Store reduced dataset (vastly smaller footprint)
4. Later: `expand_imputed_data()` - reconstruct full dataset when needed for analysis

**State Management:**

- **Immutable data** - functions do not modify input datasets, return new objects
- **Metadata preservation** - IMPID column tracks imputation identity through workflow
- **Analysis object** - custom class wrapping results list plus metadata (method, function name, delta)
- **Type preservation** - factors/character types maintained through transformations

## Key Abstractions

**Analysis Object (`class: analysis`):**
- Purpose: Encapsulates results from applying analysis function to each imputation
- Examples: Created by `analyse_mi_data()`, consumed by `rbmi::pool()`
- Pattern: S3 class with `results` (list of per-imputation results), `method`, `delta`, `fun` metadata
- Structure defined in: `R/analyse_mi_data.R` lines 264-300 (as_analysis2 function)

**Vars Object (`class: vars`):**
- Purpose: Encapsulates variable metadata for analysis (subject ID, visit, group, outcome, covariates)
- Examples: Created by `rbmi::set_vars()`
- Pattern: Named list with required fields: subjid, visit, group, outcome; optional: covariates, strategy
- Used throughout to enable flexible variable specification

**Imputation ID Column (IMPID):**
- Purpose: Identifies which imputation iteration a row belongs to
- Pattern: Numeric or character column added by `get_imputed_data()`, `create_impid()`, `expand_imputed_data()`
- Critical for `analyse_mi_data()` which splits data by IMPID and applies analysis to each subset

**Tidy Results Tibble:**
- Purpose: Publication-ready results with one row per parameter (treatment effect or LSM)
- Examples: Output from `tidy_pool_obj()`
- Pattern: Columns: parameter, description, visit, parameter_type, lsm_type, est, se, lci, uci, pval
- Structure defined in: `R/tidiers.R` lines 88-167 (tidy_pool_obj function)

## Entry Points

**validate_data():**
- Location: `R/data_helpers.R` line 64
- Triggers: Called by user before imputation
- Responsibilities: Pre-flight validation of data structure, variable existence, types, missing patterns, ICE data

**analyse_mi_data():**
- Location: `R/analyse_mi_data.R` line 79
- Triggers: Called after imputation with IMPID-tagged datasets
- Responsibilities: Apply analysis function to each imputation, construct analysis object, validate inputs

**tidy_pool_obj():**
- Location: `R/tidiers.R` line 88
- Triggers: Called on rbmi::pool() output
- Responsibilities: Transform pooled object to tidy tibble, extract visit info, add descriptions

**prepare_data_ice():**
- Location: `R/data_helpers.R` line 279
- Triggers: Called when handling intercurrent events
- Responsibilities: Convert ICE flag column to proper data_ice format for rbmi::draws()

**get_imputed_data():**
- Location: `R/utils.R` line 52
- Triggers: Called after rbmi imputation to extract datasets
- Responsibilities: Extract imputed datasets, map original subject IDs, add IMPID column

## Error Handling

**Strategy:** Fail early with detailed, actionable messages

**Patterns:**

**1. Input Validation (assertthat package):**
```R
assertthat::assert_that(
  is.data.frame(data),
  msg = "`data` must be a data.frame"
)
```
Used in: `R/analyse_mi_data.R`, `R/data_helpers.R`, `R/imputation_storage.R`, `R/result_helpers.R`

**2. Collected Issue Reporting:**
```R
issues <- character(0)
# ... check and accumulate issues ...
if (length(issues) > 0) {
  stop(paste(c("Data validation failed:", issues), collapse = "\n  - "), call. = FALSE)
}
```
Used in: `R/data_helpers.R` validate_data() - collects all problems before reporting once

**3. Type-Specific Checks with Helpful Messages:**
Example from `R/data_helpers.R`:
- Factor validation: checks class, warns on character (will be auto-converted by rbmi)
- Numeric outcome validation: explicit type check
- Covariate completeness: checks for any missing values in covariates
- Subject-visit uniqueness: ensures one row per subject-visit combination

**4. Conditional Warnings for Non-Critical Issues:**
```R
warning(sprintf("Data contains %d imputations but method expects %d", n_impids, n_expected), call. = FALSE)
```
Used in: `R/analyse_mi_data.R` lines 164-170 - when imputation count mismatch is auto-fixable

## Cross-Cutting Concerns

**Logging:** Console output via cat() and print() for S3 methods; no persistent logging framework

**Validation:** Heavy use at entry points; cascade validation pattern (data → vars → method → delta)

**Authentication:** Not applicable - analysis package, no external auth

**Variable Naming:** Uses factors and explicit variable references through vars object; avoids positional argument confusion

**Missing Data Handling:** Explicitly validates complete covariates; outcome missingness is expected/required for imputation

**Formula Construction:** Helper functions in `R/analysis_utils.R` (`as_simple_formula2()`, `extract_covariates2()`) construct formulas from variable specs

---

*Architecture analysis: 2026-02-07*
