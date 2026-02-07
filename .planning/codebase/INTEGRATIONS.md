# External Integrations

**Analysis Date:** 2026-02-07

## APIs & External Services

**Statistical Computing:**
- rbmi package - Provides reference-based multiple imputation backend
  - SDK/Client: R package (imported via `Imports: rbmi (>= 1.4)`)
  - Integration: Core; rbmiUtils wraps and extends rbmi workflows
  - Functions used: `rbmi::draws()`, `rbmi::ancova()`, `rbmi::pool()`, `rbmi::set_vars()`, `rbmi::method_bayes()`, `rbmi::method_approxbayes()`, `rbmi::method_condmean()`, `rbmi::control_bayes()`

**Binary Outcome Analysis:**
- beeca package - Marginal effects estimation for G-computation
  - SDK/Client: R package (imported via `Imports: beeca`)
  - Integration: Used in `gcomp_responder()` and `gcomp_responder_multi()` functions
  - Functions used: `beeca::get_marginal_effect()`

## Data Storage

**Databases:**
- None - Package is in-memory only; no persistent database integration

**File Storage:**
- Local filesystem only
- Data inputs: CSV format (loaded via readr in data-raw scripts)
- Data outputs: R data frames, tibbles, or user-written files

**Caching:**
- None - No caching layer; processing happens in-memory per session

## Authentication & Identity

**Auth Provider:**
- None required - Pure statistical computing; no user authentication

## Monitoring & Observability

**Error Tracking:**
- None - Package relies on standard R error handling and base warnings

**Logs:**
- Console output only via `base::message()`, `base::warning()`, and `stop()`
- No structured logging framework

## CI/CD & Deployment

**Hosting:**
- GitHub (repository: `openpharma/rbmiUtils`)
- CRAN (published package; distributed via `install.packages("rbmiUtils")`)
- GitHub Pages (documentation website at `https://openpharma.github.io/rbmiUtils/`)

**CI Pipeline:**
- GitHub Actions
- Workflows:
  - `R-CMD-check.yaml`: Package validation across platforms
  - `test-coverage.yaml`: Code coverage tracking
  - `pkgdown.yaml`: Documentation site building and deployment
  - `rhub.yaml`: R Hub cross-platform validation

**Deployment Triggers:**
- Push to main branch (R-CMD-check, pkgdown)
- Pull requests (R-CMD-check, test-coverage)
- GitHub releases (pkgdown site update)
- Manual workflow dispatch (pkgdown)

## Environment Configuration

**Required env vars:**
- None for package functionality
- GitHub Actions environment:
  - `GITHUB_PAT`: GitHub token (auto-injected, used for dependency installation)
  - `R_KEEP_PKG_SOURCE`: Set to "yes" to preserve package sources in CI

**Secrets location:**
- GitHub repository secrets: `GITHUB_PAT` (auto-available in Actions)
- No custom secrets required for package operation

## Webhooks & Callbacks

**Incoming:**
- None - Package is a library; does not expose API endpoints

**Outgoing:**
- GitHub Pages deployment hook (triggered after documentation build in `pkgdown.yaml`)
- No external webhooks configured

## Package Distribution

**CRAN:**
- Package is published on CRAN
- Installation: `install.packages("rbmiUtils")`
- Version management in `DESCRIPTION` file

**GitHub:**
- Development version available: `remotes::install_github("openpharma/rbmiUtils")`
- Releases tagged in git

## Documentation & Examples

**Data Sources:**
- Built-in datasets:
  - `ADMI`: 100,000 rows; simulated multiple imputation trial data with 12 columns (USUBJID, STRATA, REGION, TRT, BASE, CHG, AVISIT, IMPID, CRIT1FLN, CRIT1FL, CRIT)
  - `ADEFF`: Simulated efficacy dataset with missing data (loaded via CSV in data-raw)
  - Source CSV files: `data-raw/extdata/ADMI.csv` and `data-raw/extdata/ADEFF.csv`

**Vignettes:**
- `analyse2.Rmd`: Analysis workflow examples
- `data-preparation.Rmd`: Data preparation guide
- `efficient-storage.Rmd`: Data storage optimization

---

*Integration audit: 2026-02-07*
