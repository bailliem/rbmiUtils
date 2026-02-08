# Convert Pool Object to ARD Format

Converts an rbmi pool object to the pharmaverse Analysis Results Dataset
(ARD) standard using the cards package. The ARD format is a long-format
data frame where each row represents a single statistic for a given
parameter, with grouping columns for visit, parameter type, and
least-squares-mean type.

## Usage

``` r
pool_to_ard(pool_obj, conf.level = NULL)
```

## Arguments

- pool_obj:

  A pooled analysis object of class `"pool"`, typically obtained from
  [`rbmi::pool()`](https://openpharma.github.io/rbmi/latest-tag/reference/pool.html)
  after calling
  [`analyse_mi_data()`](https://openpharma.github.io/rbmiUtils/reference/analyse_mi_data.md).

- conf.level:

  Confidence level used for CI labels (e.g., `0.95` produces "95% CI
  Lower"). If `NULL` (the default), the value is taken from
  `pool_obj$conf.level`. If that is also `NULL`, defaults to `0.95`.

## Value

A data frame of class `"card"` (ARD format) with grouping columns for
visit (`group1`), parameter_type (`group2`), and lsm_type (`group3`).
Each parameter produces rows for five statistics: estimate, std.error,
conf.low, conf.high, and p.value, plus a method row.

## Details

The function works by:

1.  Tidying the pool object via
    [`tidy_pool_obj()`](https://openpharma.github.io/rbmiUtils/reference/tidy_pool_obj.md)

2.  Reshaping each parameter into long-format ARD rows (one row per
    statistic)

3.  Adding grouping columns (visit, parameter_type, lsm_type)

4.  Applying
    [`cards::as_card()`](https://insightsengineering.github.io/cards/latest-tag/reference/as_card.html)
    and
    [`cards::tidy_ard_column_order()`](https://insightsengineering.github.io/cards/latest-tag/reference/tidy_ard_order.html)
    for standard ARD structure

The resulting ARD passes
[`cards::check_ard_structure()`](https://insightsengineering.github.io/cards/latest-tag/reference/check_ard_structure.html)
validation and is suitable for downstream use with gtsummary.

## See also

[`tidy_pool_obj()`](https://openpharma.github.io/rbmiUtils/reference/tidy_pool_obj.md),
[`cards::as_card()`](https://insightsengineering.github.io/cards/latest-tag/reference/as_card.html),
[`cards::check_ard_structure()`](https://insightsengineering.github.io/cards/latest-tag/reference/check_ard_structure.html)

## Examples

``` r
# \donttest{
# Requires the cards package
if (requireNamespace("cards", quietly = TRUE)) {
  # After running an rbmi analysis pipeline:
  # pool_obj <- rbmi::pool(analysis_obj)
  # ard <- pool_to_ard(pool_obj)
  # cards::check_ard_structure(ard)
}
#> NULL
# }
```
