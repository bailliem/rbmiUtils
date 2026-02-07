# Codebase Concerns

**Analysis Date:** 2026-02-07

## Test Coverage Gaps

**Internal functions not testable by default:**
- Files: `R/analysis_utils.R`
- What's not tested: Helper functions `extract_covariates2()` and `as_simple_formula2()` are marked as internal (`@keywords internal`) but test suite attempts to call them directly, causing test failures
- Files affected: `tests/testthat/test-analysis_utils.R` lines 2-20
- Risk: Test suite failures (currently 8 tests error out) reduce confidence in code quality. Internal functions bypass validation in wrapper functions.
- Priority: Medium
- Recommendation: Either export and document these utilities, or refactor tests to only test through public API (`gcomp_responder()`). Consider whether these utilities should be public for user extension.

**Data package loading in tests:**
- Files: All test files reference `ADMI` and `ADEFF` datasets
- What's not tested: Multiple tests fail with "data set 'ADMI' not found" warnings, indicating tests cannot load package datasets during automated testing
- Test files: `tests/testthat/test-*.R` (affects at least 14 tests)
- Risk: Tests cannot run without manual data loading; CI/CD pipelines may fail or skip tests. Reduces reliability of automated testing.
- Recommendation: Add `data(ADMI); data(ADEFF)` within test functions, or use `testthat::skip_if_not()` to handle missing data gracefully.

**G-computation test failures:**
- Files: `R/analysis_utils.R`
- Functions affected: `gcomp_responder()`, `gcomp_responder_multi()`
- Issue: Multiple tests attempt to call non-exported utility functions that don't exist in the namespace
- Risk: Test suite cannot validate the core g-computation workflow for binary outcomes
- Priority: High
- Fix approach: Export utility functions or refactor internal dependencies

## Type Coercion and Edge Cases

**String conversion in matching logic:**
- Files: `R/imputation_storage.R` lines 120-134
- Issue: Heavy reliance on `as.character()` conversions for robust matching with separator `"|||\`". While robust, this approach could mask type mismatches that should be errors (e.g., matching numeric factors as strings).
- Risk: Silent failures if data types change unexpectedly; reduces early error detection
- Impact: `reduce_imputed_data()` may silently produce wrong results if key columns change type
- Improvement path: Add explicit type checks before conversion; document expected column types in error messages

**Factor level handling inconsistency:**
- Files: `R/data_helpers.R` lines 369-372
- Issue: `prepare_data_ice()` treats factors and non-factors differently when determining visit order. Non-factor visits fall back to sequence order rather than preserving any meaningful order.
- Risk: Intermittent bugs if visit columns aren't consistently factored. Users may miss visits if data isn't properly ordered.
- Safe modification: Enforce factor requirement in validation; document that visit column must be ordered factor

**Character/Factor coercion in matching:**
- Files: `R/data_helpers.R` lines 488-496
- Issue: In `summarise_missingness()`, subject and visit comparisons repeatedly convert to character with `as.character()`. For large datasets, this creates many temporary strings.
- Risk: Performance degradation with large N subjects or visits. Memory inefficiency.
- Improvement: Use integer indices or factors directly; only convert for comparison when necessary

## Performance Bottlenecks

**Repeated string conversions in loops:**
- Files: `R/imputation_storage.R` line 131-132, `R/expand_imputed_data.R` lines 266-269
- Issue: String concatenation with `paste()` using `"|||\`" separator is called on every row during matching operations
- Risk: O(n) string allocations in already-intensive operations. Performance impact grows with dataset size.
- Current scope: `reduce_imputed_data()` and `expand_imputed_data()` with datasets containing thousands of imputations
- Scaling limit: Likely bottleneck when working with >100k rows or very large number of imputations
- Improvement: Use `interaction()` or merge on index columns instead of string concatenation

**Lapply with data frame subsetting:**
- Files: `R/analyse_mi_data.R` lines 218-226
- Issue: `data |> group_split(IMPID) |> lapply()` creates separate copies of large data frames for each imputation. No vectorization.
- Risk: Memory usage scales linearly with number of imputations × dataset size. For 1000 imputations of large datasets, can exhaust RAM.
- Current usage: Acceptable for typical workflow (10-100 imputations), breaks at >500 imputations
- Recommendation: Document this limitation; provide guidance on maximum safe imputation counts

**Repeated attribute restoration:**
- Files: `R/imputation_storage.R` lines 138-146 and `R/expand_imputed_data.R` lines 313-324
- Issue: Column attributes (e.g., labels, formats) are restored in nested loops, rewriting same attributes repeatedly
- Risk: Unnecessary overhead. Notable with datasets having 50+ columns.
- Impact: Measurable slowdown in `reduce_imputed_data()` and `expand_imputed_data()` with wide datasets
- Improvement: Collect attributes before loop; apply once at end

## Data Validation Gaps

**No validation of subject-visit uniqueness in delta data:**
- Files: `R/analyse_mi_data.R` lines 195-214
- Issue: `delta` parameter validated for column presence and variable names, but not checked for duplicate subject-visit combinations
- Risk: Silent failure if delta data has duplicates; merge will replicate rows or produce unexpected results
- Safe modification: Add check: `if (nrow(delta) != nrow(unique(delta[, c(vars$subjid, vars$visit)]))) stop(...)`

**Character column handling in prepare_data_ice():**
- Files: `R/data_helpers.R` lines 308-312
- Issue: If `vars$strategy` is NULL or empty, function silently defaults to `strategy = "strategy"` (string literal). No warning to user.
- Risk: User may not realize their strategy column is being created with wrong name, leading to downstream errors
- Priority: Medium
- Recommendation: Throw error instead of silently defaulting; require explicit strategy specification

**Covariate interaction term parsing:**
- Files: `R/data_helpers.R` lines 92-97
- Issue: Interaction terms in covariates are split naively with regex `"[\\*\\:]"` without validation that all split components are real columns
- Risk: If interaction term is malformed (e.g., `"A*"` or `"A:"`), parsing succeeds but later steps fail with confusing error about missing column
- Impact: Hard-to-debug user errors when specifying covariates
- Recommendation: Validate split terms are non-empty after trim; provide early error with clear guidance

## Fragile Areas

**tidy_pool_obj() parameter parsing:**
- Files: `R/tidiers.R` lines 98-106
- What makes it fragile: Function splits `parameter` column on `_` separator without validation. If analysis functions produce parameter names with underscores (e.g., `"trt_Drug_vs_Placebo_Week24"`), separator becomes ambiguous.
- Safe modification: Document naming convention strictly; validate exactly 3 parts expected; throw error if splitting produces unexpected count
- Test coverage: Parameter splitting logic not explicitly tested; only indirect testing through integration tests

**Attribute preservation across rbmi pipeline:**
- Files: Multiple files (`imputation_storage.R`, `expand_imputed_data.R`)
- Fragility: Attributes (labels, classes, formats) are selectively restored using hardcoded exclude lists (`setdiff(names(col_attrs), c("names", "dim", "dimnames"))`)
- Risk: If rbmi or downstream tools add new standard attributes, they may be lost or incorrectly handled
- Impact: Loss of metadata (e.g., variable labels, units) during imputation
- Recommendation: Document which attributes are preserved; provide control via parameter

**beeca integration in gcomp_responder():**
- Files: `R/analysis_utils.R` lines 72-81
- Fragility: Tightly coupled to beeca package's output format and column names (`STAT`, `STATVAL`, `TRTVAL`)
- Risk: Breaking changes in beeca versions could fail silently or produce cryptic errors
- Mitigation present: None (no version pinning in DESCRIPTION)
- Recommendation: Add version constraint on beeca dependency (currently >= 1.4 on rbmi, but no beeca version specified)

## Dependencies at Risk

**Underdocumented dependency on beeca output format:**
- Package: beeca
- Files: `R/analysis_utils.R` lines 83-100
- Risk: Code directly accesses `beeca::get_marginal_effect()` output structure (columns: STAT, STATVAL, TRTVAL). No validation of output format.
- Impact: Upgrade to new beeca version could break g-computation for binary outcomes without warning
- Recommendation: Add version constraint to DESCRIPTION; add unit tests that validate beeca output schema

**rbmi version compatibility:**
- Package: rbmi
- Constraint: `>= 1.4`
- Risk: Code in `as_analysis2()` switches on `class(method)[[2]]` to detect method type. If rbmi changes class hierarchy, breaks silently.
- Files: `R/analyse_mi_data.R` lines 148-155, 271-281
- Recommendation: Add explicit version pinning or add validation that method classes match expected types

**Dependency on tidyr separate behavior:**
- Package: tidyr
- Files: `R/tidiers.R` lines 98-106
- Issue: `tidyr::separate()` with `fill = "right"` and `extra = "merge"` behavior may change in future versions
- Risk: Parameter splitting in `tidy_pool_obj()` could behave unexpectedly with future tidyr
- Mitigation: Moderate (tidyr stable), but undocumented fragility

## Security Considerations

**Formula construction from user input:**
- Files: `R/analysis_utils.R` lines 205-214
- Risk: `as_simple_formula2()` builds formulas from user-supplied covariate strings without validation
- Attack vector: User provides formula injection via interaction terms
- Current mitigation: assertthat checks in wrapper functions, but this internal utility has none
- Recommendation: Add input sanitization for covariate names (allow only alphanumeric, underscore, dot)

## Missing Critical Features

**No explicit imputation method validation:**
- Files: `R/analyse_mi_data.R` lines 132-134
- Issue: Function checks if method is NULL but doesn't validate that it's a recognized rbmi method object
- Risk: User passes invalid method object; error occurs deep in analysis function with unclear origin
- Recommendation: Add explicit type check: `stopifnot(inherits(method, "method_bayes") | inherits(method, "method_condmean") | ...)`

**No warning for truncated imputations:**
- Files: `R/analyse_mi_data.R` lines 161-173
- Issue: When data contains more imputations than expected by method, warning is issued but data is silently filtered. No record of which imputations were discarded.
- Risk: Silent data loss; user may not realize analyses are based on subset of data
- Recommendation: Also return a list of discarded IMPID values for user inspection

**No support for custom pooling methods:**
- Scope: `analyse_mi_data()` -> `as_analysis2()` -> hardcoded pooling class assignment
- Files: `R/analyse_mi_data.R` lines 271-281
- Risk: Users cannot add custom pooling methods without modifying package code
- Recommendation: Consider allowing `pooling_method` parameter to override default assignment

## Known Warnings

**Character columns converted to factors warning:**
- Source: `R/data_helpers.R` line 120-123
- Issue: `validate_data()` warns if character columns will be converted to factors by `rbmi::draws()`
- Risk: Warnings can hide real issues if user expects silent coercion
- Mitigation: Warning clearly explains behavior; users can convert explicitly before validation

**R-hat convergence warnings in tests:**
- Source: Test output shows "largest R-hat is 1.06, indicating chains have not mixed"
- Issue: Bayesian imputation tests use limited warmup (20-200 iterations) for speed; produces convergence warnings
- Risk: New developers may misinterpret as code failure rather than test configuration
- Recommendation: Document this is expected in quick tests; use longer chains for production validation

## Deprecation Warnings

**Deprecated `stringsAsFactors` parameter:**
- Files: `R/data_helpers.R` lines 359, 387, and multiple locations
- Issue: `stringsAsFactors = FALSE` in data.frame() calls is deprecated in R 4.1+, though still functional
- Impact: Code future-proofs itself against R 5.0 when parameter removal is likely
- Recommendation: In next major release, remove this parameter entirely; rely on R defaults

## Technical Debt Summary

| Area | Severity | Impact | Files | Fix Effort |
|------|----------|--------|-------|-----------|
| Test failures on internal functions | High | Cannot run test suite | `test-analysis_utils.R` | Medium |
| String-based matching performance | Medium | Slow with large data | `imputation_storage.R` | High |
| Memory usage in lapply loop | Medium | Crashes with many imputations | `analyse_mi_data.R` | Medium |
| Parameter parsing fragility | Medium | Silent failures | `tidiers.R` | Low |
| Missing delta validation | Low | Silent wrong results | `analyse_mi_data.R` | Low |
| Undocumented beeca dependency | Medium | Version incompatibility | `analysis_utils.R` | Medium |

---

*Concerns audit: 2026-02-07*
