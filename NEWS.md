# rbmiUtils 0.2.0

## New Features

* End-to-end pipeline vignette "From rbmi Analysis to Regulatory Tables" with
  continuous and binary analysis walkthroughs.
* README enhanced with rendered efficacy table and forest plot visual teasers.
* Rendered example output images in `plot_forest()` and `efficacy_table()` help pages.
* Inline cross-references to rbmi and beeca documentation in all vignettes.

## Improvements

* `validate_data()` now uses cli-formatted error messages with clearer guidance
  for malformed interaction terms, empty data, and type mismatches.
* `prepare_data_ice()` now errors immediately when `vars$strategy` is NULL
  instead of silently using a default column name.
* `prepare_data_ice()` warns when visit column is character with guidance to
  convert to factor for correct ordering.
* `validate_data()` batches all type coercion warnings into a single message
  and warns on all-NA covariate columns.
* Data preparation functions handle edge cases (single subject, single visit,
  all-NA outcome, all-complete data) gracefully.

# rbmiUtils 0.1.0

## Breaking Changes

* `tidy_pool_obj()` now uses regex-based parameter parsing instead of splitting
  on `_`. Output columns (`parameter_type`, `lsm_type`, `visit`) now contain
  the full visit name rather than truncated fragments.

## New Features

* `efficacy_table()` creates regulatory-style gt summary tables from pool objects.
* `plot_forest()` creates publication-quality three-panel forest plots from pool objects.
* `pool_to_ard()` converts pool objects to pharmaverse ARD format via the cards package.
* `print()` and `summary()` S3 methods for pool and analysis objects.
* `create_impid()` converts lists of imputed datasets into stacked data.frames.
* `combine_results()` combines tidy results from multiple analyses.
* `format_results()` and `format_results_table()` for publication-ready formatting.
* `extract_trt_effects()` and `extract_lsm()` convenience filters for tidy results.

## Dependencies

* Added `cli` (>= 3.6.0) and `lifecycle` (>= 1.0.4) to Imports for improved
  error messaging and deprecation support.

## Improvements

* Enhanced input validation in `analyse_mi_data()` with clearer error messages
  using `inherits()` instead of internal class helpers.
* Enhanced input validation in `prepare_data_ice()` to check vars fields.
* Cross-references between related functions in documentation.
* Expanded test coverage and integration tests.

## Previous Releases

* `reduce_imputed_data()` and `expand_imputed_data()` for efficient storage
  of imputed datasets.
* `validate_data()` for pre-flight validation before `rbmi::draws()`.
* `prepare_data_ice()` to build `data_ice` from flagged ICE columns.
* `summarise_missingness()` for tabulating missing data patterns.
* `format_pvalue()`, `format_estimate()`, `format_results_table()` formatting utilities.
* `ADEFF` and `ADMI` example datasets.
* Initial documentation via pkgdown.
