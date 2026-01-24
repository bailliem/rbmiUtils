# G-computation Analysis for a Single Visit

Performs logistic regression and estimates marginal effects for binary
outcomes.

## Usage

``` r
gcomp_responder(
  data,
  vars,
  reference_levels = NULL,
  var_method = "Ge",
  type = "HC0",
  contrast = "diff"
)
```

## Arguments

- data:

  A data.frame with one visit of data.

- vars:

  A list containing `group`, `outcome`, `covariates`, and `visit`.

- reference_levels:

  Optional vector specifying reference level(s) of the treatment factor.

- var_method:

  Marginal variance estimation method (default: "Ge").

- type:

  Type of robust variance estimator (default: "HC0").

- contrast:

  Type of contrast to compute (default: "diff").

## Value

A named list containing estimates and standard errors for treatment
comparisons and within-arm means.
