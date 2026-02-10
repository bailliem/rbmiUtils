# rbmiUtils

`rbmiUtils` bridges [rbmi](https://github.com/openpharma/rbmi) analysis
results into publication-ready regulatory tables and forest plots. It
extends rbmi for clinical trial workflows, handling everything from data
validation through to formatted efficacy outputs.

## Installation

You can install the package from CRAN or the development version from
GitHub:

| Type        | Source | Command                                           |
|-------------|--------|---------------------------------------------------|
| Release     | CRAN   | `install.packages("rbmiUtils")`                   |
| Development | GitHub | `remotes::install_github("openpharma/rbmiUtils")` |

## Quick Start

Starting from an [rbmi](https://insightsengineering.github.io/rbmi/)
pool object, `rbmiUtils` produces publication-ready outputs in a few
lines:

``` r
library(rbmiUtils)
library(rbmi)

# Analyse across imputations and pool results
ana_obj  <- analyse_mi_data(data = ADMI, vars = vars, method = method)
pool_obj <- pool(ana_obj)

# Publication outputs
efficacy_table(pool_obj, arm_labels = c(ref = "Placebo", alt = "Drug A"))
plot_forest(pool_obj, arm_labels = c(ref = "Placebo", alt = "Drug A"))
```

**Forest Plot**

![Forest Plot](reference/figures/README-forest-plot-1.png)

Forest Plot

**Efficacy Table**

![Efficacy Table](reference/figures/README-efficacy-table-1.png)

Efficacy Table

See the [end-to-end pipeline
vignette](https://openpharma.github.io/rbmiUtils/articles/pipeline.html)
for the complete walkthrough from raw data to these outputs.

## Key Features

- [`validate_data()`](https://openpharma.github.io/rbmiUtils/reference/validate_data.md)
  – pre-flight checks on data structure before imputation
- [`analyse_mi_data()`](https://openpharma.github.io/rbmiUtils/reference/analyse_mi_data.md)
  – run ANCOVA (or custom analysis) across all imputations
- [`tidy_pool_obj()`](https://openpharma.github.io/rbmiUtils/reference/tidy_pool_obj.md)
  – tidy pooled results with visit-level annotations
- [`efficacy_table()`](https://openpharma.github.io/rbmiUtils/reference/efficacy_table.md)
  – regulatory-style gt tables (CDISC/ICH Table 14.2.x format)
- [`plot_forest()`](https://openpharma.github.io/rbmiUtils/reference/plot_forest.md)
  – three-panel forest plots with estimates, CIs, and p-values
- [`pool_to_ard()`](https://openpharma.github.io/rbmiUtils/reference/pool_to_ard.md)
  – convert pool objects to pharmaverse ARD format
- [`get_imputed_data()`](https://openpharma.github.io/rbmiUtils/reference/get_imputed_data.md)
  – extract long-format imputed datasets
- [`format_pvalue()`](https://openpharma.github.io/rbmiUtils/reference/format_pvalue.md)
  /
  [`format_estimate()`](https://openpharma.github.io/rbmiUtils/reference/format_estimate.md)
  – publication-ready formatting

## Learn More

- [From rbmi Analysis to Regulatory
  Tables](https://openpharma.github.io/rbmiUtils/articles/pipeline.html)
  – end-to-end walkthrough from raw data to regulatory outputs
- [Storing and Analyzing Imputed
  Data](https://openpharma.github.io/rbmiUtils/articles/analyse2.html) –
  focused guide on analysis workflows
- [Package documentation](https://openpharma.github.io/rbmiUtils/)

## Development Status

This package is experimental and under active development. Feedback and
contributions are welcome via [GitHub
issues](https://github.com/openpharma/rbmiUtils/issues) or pull
requests.
