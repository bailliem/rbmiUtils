## R CMD check results

0 errors | 0 warnings | 0 notes

## Summary of Submission (v0.3.0)

This is a re-submission with minor improvements 

* `describe_draws()` extracts structured metadata from rbmi draws objects,
  including method type, formula, sample count, and (for Bayesian methods)
  MCMC convergence diagnostics (ESS, Rhat).
* `describe_imputation()` extracts imputation metadata including method,
  number of imputations (M), reference arm mappings, and a missingness
  breakdown by visit and treatment arm.
* `pool_to_ard()` gains an `analysis_obj` parameter that enriches the ARD
  with MI diagnostic statistics (FMI, lambda, RIV, Barnard-Rubin adjusted
  df, relative efficiency) when the pooling method is Rubin's rules.

## Test environments

The package was tested in the following environments:

- macOS, R release (Local Machine)
- Windows, R release (Local Machine)
- Fedora, R release (Local Machine)
- macOS, devel (macOS builder)
