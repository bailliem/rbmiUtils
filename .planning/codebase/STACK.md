# Technology Stack

**Analysis Date:** 2026-02-07

## Languages

**Primary:**
- R (>= 4.1) - All package code, functions, and examples

## Runtime

**Environment:**
- R >= 4.1.0 (specified in `DESCRIPTION`)

**Package Manager:**
- R package system (DESCRIPTION for dependencies)
- Lockfile: Not applicable (R packages use DESCRIPTION/NAMESPACE)

## Frameworks

**Core:**
- rbmi (>= 1.4) - Reference-based multiple imputation; foundational dependency providing imputation, pooling, and analysis infrastructure

**Data Manipulation:**
- dplyr - Data frame manipulation and pipelines
- purrr - Functional programming utilities (map, flatten operations)
- tidyr - Data tidying (separate, pivot operations)
- rlang - Quasiquotation and NSE support (`.data` pronoun usage)

**Validation:**
- assertthat - Pre-condition assertions and input validation

**Testing:**
- testthat (>= 3.0.0) - Unit testing framework (Config: `Config/testthat/edition: 3`)

**Documentation:**
- roxygen2 (version 7.3.2) - Source code documentation from Roxygen comments (auto-generates NAMESPACE and .Rd files)
- knitr - Vignette rendering
- rmarkdown - R Markdown support for README.Rmd and vignettes

**Statistical/Specialized:**
- beeca - Marginal effects estimation for binary outcomes (G-computation); used in `gcomp_responder()` and `gcomp_responder_multi()`
- rstan - Optional: Bayesian modeling support (suggested, used indirectly via rbmi for Bayesian imputation)

**Supplementary (Suggested, not required):**
- readr - CSV reading (used only in data-raw scripts)
- tibble - Tibble display and manipulation (used in examples and output formatting)
- ggplot2 - Visualization (optional for user workflows)
- spelling - Spell checking (optional development tool)

## Build & Documentation

**Documentation:**
- Configuration file: `DESCRIPTION` (metadata, dependencies)
- Roxygen control: `Roxygen: list(markdown = TRUE)` enables markdown in roxygen comments
- Build system: `PackageRoxygenize: rd,collate,namespace` (auto-generates .Rd, NAMESPACE, collation order)

**Package Site:**
- pkgdown - Website generation
- Configuration: `_pkgdown.yml` with bootstrap v5 template
- Deployment: GitHub Pages (via GitHub Actions)

**R Project Settings:**
- Configuration file: `rbmiUtils.Rproj`
- Roxygen-based code generation
- Package build: `--no-multiarch --with-keep.source`

## Key Dependencies

**Critical:**
- **rbmi** (>= 1.4) - Core imputation and pooling framework; without this the package cannot function
- **dplyr** - Data manipulation pipeline; used throughout for filtering, selecting, grouping, and piping
- **beeca** - G-computation for binary outcomes; enables responder analysis with marginal effects estimation
- **assertthat** - Input validation; prevents silent failures from malformed data or parameters

**Infrastructure:**
- **tidyr** - Separating and restructuring parameter columns for tidying pooled results
- **purrr** - Functional iteration (group_map, flatten) for result processing
- **rlang** - NSE support and data masking (`.data` pronoun in dplyr operations)

## Configuration

**Environment:**
- No environment variables required for base functionality
- GitHub Actions use `GITHUB_PAT` for accessing private repositories (injected at runtime)

**Build Configuration:**
- Ignored from distribution: `.Rbuildignore` excludes:
  - R project files (`.Rproj.user`)
  - README source (`.Rmd`)
  - Data generation scripts (`data-raw`)
  - Documentation outputs (`docs`, `Meta`)
  - GitHub config (`.github`, `.nojekyll`)
  - pkgdown config (`_pkgdown.yml`)

## Platform Requirements

**Development:**
- R >= 4.1
- Roxygen2 (for documentation generation)
- Pandoc (via r-lib/actions in CI)
- RcppEigen, BH (for rstan dependency support in CI matrix)

**Testing:**
- Multi-platform CI matrix:
  - macOS latest (R release)
  - Windows latest (R release)
  - Ubuntu latest (R release, devel, oldrel-1)

**Production/Deployment:**
- GitHub Pages (static site hosting for package documentation)
- CRAN distribution (package is published on CRAN)
- No server runtime required; pure R package distribution

## CI/CD

**Build Verification:**
- Workflow: `.github/workflows/R-CMD-check.yaml`
- Runs R CMD check across multiple platforms and versions
- Extra packages installed: rcmdcheck, BH, RcppEigen, rstan, QuickJSR
- Build args: `--no-manual --compact-vignettes=gs+qpdf`

**Test Coverage:**
- Workflow: `.github/workflows/test-coverage.yaml`
- Tracks code coverage metrics

**Documentation:**
- Workflow: `.github/workflows/pkgdown.yaml`
- Builds and deploys package website to GitHub Pages
- Deployed on: main branch push, release publication, manual dispatch

**Cross-Platform Validation:**
- Workflow: `.github/workflows/rhub.yaml`
- Validates against R Hub service for additional platform checks

---

*Stack analysis: 2026-02-07*
