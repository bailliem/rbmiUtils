# G-computation for a Binary Outcome at Multiple Visits

Applies
[`gcomp_responder()`](https://openpharma.github.io/rbmiUtils/reference/gcomp_responder.md)
separately for each unique visit in the data.

## Usage

``` r
gcomp_responder_multi(data, vars, reference_levels = NULL, ...)
```

## Arguments

- data:

  A data.frame containing multiple visits.

- vars:

  A list specifying analysis variables.

- reference_levels:

  Optional reference level for the treatment variable.

- ...:

  Additional arguments passed to
  [`gcomp_responder()`](https://openpharma.github.io/rbmiUtils/reference/gcomp_responder.md).

## Value

A named list of estimates for each visit and treatment group.

## Examples

``` r
# \donttest{
library(dplyr)
library(rbmi)
library(rbmiUtils)

data("ADMI")

ADMI <- ADMI |>
  mutate(
    TRT = factor(TRT, levels = c("Placebo", "Drug A")),
    STRATA = factor(STRATA),
    REGION = factor(REGION)
  )

# Note: method must match the original used for imputation
method <- method_bayes(
  n_samples = 100,
  control = control_bayes(warmup = 20, thin = 2)
)

vars_binary <- set_vars(
  subjid = "USUBJID",
  visit = "AVISIT",
  group = "TRT",
  outcome = "CRIT1FLN",
  covariates = c("BASE", "STRATA", "REGION")
)

ana_obj_prop <- analyse_mi_data(
  data = ADMI,
  vars = vars_binary,
  method = method,
  fun = gcomp_responder_multi,
  reference_levels = "Placebo",
  contrast = "diff",
  var_method = "Ge",
  type = "HC0"
)

pool(ana_obj_prop)
#> 
#> ── Pool Object ─────────────────────────────────────────────────────────────────
#> 6 parameters across 4 visits
#> Method: rubin
#> N imputations: 100
#> Confidence: 95%
#> ────────────────────────────────────────────────────────────────────────────────
#>                   parameter                  visit   est   lci   uci    pval
#>  trt_Drug A-Placebo_Week 24 Drug A-Placebo_Week 24 -0.03 -0.05 -0.01   0.007
#>          lsm_Drug A_Week 24                Week 24  0.00  0.00  0.00   0.921
#>         lsm_Placebo_Week 24                Week 24  0.03  0.01  0.05   0.007
#>  trt_Drug A-Placebo_Week 48 Drug A-Placebo_Week 48 -0.10 -0.14 -0.05 < 0.001
#>          lsm_Drug A_Week 48                Week 48  0.01  0.00  0.02   0.150
#>         lsm_Placebo_Week 48                Week 48  0.10  0.06  0.14 < 0.001
# }
```
