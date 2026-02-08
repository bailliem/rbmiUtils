# Print Method for Pool Objects

Displays a formatted summary of a pooled analysis object from
[`rbmi::pool()`](https://openpharma.github.io/rbmi/latest-tag/reference/pool.html).
Uses cli formatting to show rounded estimates, confidence intervals,
parameter labels, method information, number of imputations, and
confidence level.

## Usage

``` r
# S3 method for class 'pool'
print(x, digits = 2, ...)
```

## Arguments

- x:

  An object of class `pool`, typically obtained from
  [`rbmi::pool()`](https://openpharma.github.io/rbmi/latest-tag/reference/pool.html).

- digits:

  Integer. Number of decimal places for rounding estimates, standard
  errors, and confidence interval bounds. Default is 2.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the original pool object `x` (for pipe chaining).

## Details

This method overrides
[`rbmi::print.pool()`](https://openpharma.github.io/rbmi/latest-tag/reference/pool.html)
to provide enhanced, formatted console output using the cli package. The
override produces a "Registered S3 method overwritten" message at
package load time, which is expected and harmless (same pattern as
[`print.analysis()`](https://openpharma.github.io/rbmiUtils/reference/print.analysis.md)).

The output includes:

- A header with parameter and visit counts

- Metadata: pooling method, number of imputations, confidence level

- A compact results table with key columns: parameter, visit, est, lci,
  uci, pval

## See also

- [`tidy_pool_obj()`](https://openpharma.github.io/rbmiUtils/reference/tidy_pool_obj.md)
  for full tidy tibble output

- [`summary.pool()`](https://openpharma.github.io/rbmiUtils/reference/summary.pool.md)
  for visit-level breakdown with significance flags

- [`rbmi::pool()`](https://openpharma.github.io/rbmi/latest-tag/reference/pool.html)
  to create pool objects

## Examples

``` r
# \donttest{
library(rbmi)
library(rbmiUtils)
data("ADMI")

ADMI$TRT <- factor(ADMI$TRT, levels = c("Placebo", "Drug A"))
ADMI$USUBJID <- factor(ADMI$USUBJID)
ADMI$AVISIT <- factor(ADMI$AVISIT)

vars <- set_vars(
  subjid = "USUBJID", visit = "AVISIT", group = "TRT",
  outcome = "CHG", covariates = c("BASE", "STRATA", "REGION")
)
method <- method_bayes(n_samples = 20, control = control_bayes(warmup = 20))

ana_obj <- analyse_mi_data(ADMI, vars, method, fun = ancova)
#> Warning: Data contains 100 imputations but method expects 20. Using first 20
#> imputations.
pool_obj <- pool(ana_obj)
print(pool_obj)
#> 
#> ── Pool Object ─────────────────────────────────────────────────────────────────
#> 6 parameters across 2 visits
#> Method: rubin
#> N imputations: 20
#> Confidence: 95%
#> ────────────────────────────────────────────────────────────────────────────────
#>        parameter   visit   est   lci   uci    pval
#>      trt_Week 24 Week 24 -2.18 -2.54 -1.82 < 0.001
#>  lsm_ref_Week 24 Week 24  0.09 -0.17  0.34   0.514
#>  lsm_alt_Week 24 Week 24 -2.10 -2.34 -1.85 < 0.001
#>      trt_Week 48 Week 48 -3.79 -4.30 -3.29 < 0.001
#>  lsm_ref_Week 48 Week 48  0.04 -0.33  0.40   0.846
#>  lsm_alt_Week 48 Week 48 -3.76 -4.11 -3.41 < 0.001
# }
```
