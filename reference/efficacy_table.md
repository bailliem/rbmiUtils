# Create Regulatory-Style Efficacy Summary Table

Takes an rbmi pool object and produces a publication-ready gt table in
the style of CDISC/ICH Table 14.2.x. The table displays least squares
means by treatment arm, treatment differences, confidence intervals, and
p-values, organized by visit row groups.

## Usage

``` r
efficacy_table(
  pool_obj,
  title = NULL,
  subtitle = NULL,
  digits = 2,
  ci_level = NULL,
  arm_labels = NULL,
  pval_digits = 3,
  pval_threshold = 0.001,
  ...
)
```

## Arguments

- pool_obj:

  A pooled analysis object of class `"pool"`, typically obtained from
  [`rbmi::pool()`](https://openpharma.github.io/rbmi/latest-tag/reference/pool.html)
  after calling
  [`analyse_mi_data()`](https://openpharma.github.io/rbmiUtils/reference/analyse_mi_data.md).

- title:

  Optional character string for the table title.

- subtitle:

  Optional character string for the table subtitle.

- digits:

  Integer. Number of decimal places for estimates and standard errors.
  Default is 2.

- ci_level:

  Numeric. Confidence level for CI column labeling. If `NULL` (the
  default), extracted from `pool_obj$conf.level`. Falls back to 0.95 if
  neither is available.

- arm_labels:

  Named character vector with elements `"ref"` and `"alt"` providing
  custom labels for the reference and treatment arms. If `NULL` (the
  default), uses `"Reference"` and `"Treatment"`.

- pval_digits:

  Integer. Number of decimal places for p-values. Default is 3.

- pval_threshold:

  Numeric. P-values below this threshold are displayed as "\<
  threshold". Default is 0.001.

- ...:

  Additional arguments passed to
  [`gt::gt()`](https://gt.rstudio.com/reference/gt.html).

## Value

A gt table object of class `gt_tbl`.

## Details

This function assumes a single-parameter-per-visit pool object (the
standard output from an rbmi ANCOVA or MMRM pipeline). It internally
calls
[`tidy_pool_obj()`](https://openpharma.github.io/rbmiUtils/reference/tidy_pool_obj.md)
to parse the pool object, then constructs the gt table.

**Arm labels:** Use the `arm_labels` parameter to customize arm names in
the table. For example,
`arm_labels = c(ref = "Placebo", alt = "Drug A")` will display "LS Mean
(Placebo)" and "LS Mean (Drug A)" instead of the defaults.

**Customization:** The returned gt object can be further customized
using standard gt piping, e.g.,
`efficacy_table(pool_obj) |> gt::tab_options(...)`.

**Example output:**

![](figures/efficacy_table-example.png)

## See also

- [`tidy_pool_obj()`](https://openpharma.github.io/rbmiUtils/reference/tidy_pool_obj.md)
  for the underlying data transformation

- [`format_pvalue()`](https://openpharma.github.io/rbmiUtils/reference/format_pvalue.md)
  for p-value formatting rules

- [`rbmi::pool()`](https://openpharma.github.io/rbmi/latest-tag/reference/pool.html)
  to create pool objects

## Examples

``` r
# \donttest{
if (requireNamespace("gt", quietly = TRUE)) {
  # After running an rbmi analysis pipeline:
  # pool_obj <- rbmi::pool(analysis_obj)
  # tbl <- efficacy_table(pool_obj)
  # tbl  # renders in viewer
  #
  # With custom labels:
  # efficacy_table(pool_obj,
  #   title = "Table 14.2.1",
  #   subtitle = "ANCOVA of Change from Baseline",
  #   arm_labels = c(ref = "Placebo", alt = "Drug A")
  # )
}
#> NULL
# }
```
