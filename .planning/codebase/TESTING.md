# Testing Patterns

**Analysis Date:** 2026-02-07

## Test Framework

**Runner:**
- testthat (>= 3.0.0)
- Config: `tests/testthat.R` (standard setup)
- R package testing edition: 3 (Config/testthat/edition: 3 in DESCRIPTION)

**Assertion Library:**
- testthat built-in expectations (`expect_equal()`, `expect_error()`, etc.)

**Run Commands:**
```bash
# Run all tests (from R console)
devtools::test()

# Or from command line
R CMD check .

# Run specific test file
devtools::test(filter = "formatting")

# Watch mode (interactive development)
devtools::load_all()
# Then manually source test files or use testthat::test_file()
```

## Test File Organization

**Location:**
- Co-located in `tests/testthat/` directory (not beside source files)
- Separate from `R/` directory

**Naming:**
- Pattern: `test-{module}.R`
- Examples: `test-formatting.R`, `test-data_helpers.R`, `test-utils.R`

**Structure:**
```
tests/
├── testthat.R                          # Standard test runner config
└── testthat/
    ├── test-analyse_mi_data.R          # 405 lines, 20+ test blocks
    ├── test-analysis_utils.R           # 432 lines
    ├── test-data_helpers.R             # 437 lines
    ├── test-formatting.R               # 257 lines
    ├── test-imputation_storage.R       # 379 lines
    ├── test-integration.R              # 419 lines (end-to-end workflows)
    ├── test-result_helpers.R           # 256 lines
    ├── test-tidiers.R                  # 338 lines
    └── test-utils.R                    # 259 lines
```

**Total Test Coverage:**
- 3,182 lines of test code across 9 test files
- 159 `test_that()` blocks (approximately 17 per test file)

## Test Structure

**Suite Organization:**
```r
# From tests/testthat/test-formatting.R
test_that("format_pvalue handles basic cases", {
  # Standard p-values
  expect_equal(format_pvalue(0.5), "0.500")
  expect_equal(format_pvalue(0.05), "0.050")
  expect_equal(format_pvalue(0.001), "0.001")

  # Below threshold
  expect_equal(format_pvalue(0.0001), "< 0.001")
  expect_equal(format_pvalue(0.0009), "< 0.001")

  # Boundary cases
  expect_equal(format_pvalue(0), "< 0.001")
  expect_equal(format_pvalue(1), "1.000")
})
```

**Patterns:**

**Setup Pattern:**
- Fixture functions at top of test files (not setup blocks)
- Example from `test-data_helpers.R`:
```r
make_test_data <- function() {
  data.frame(
    USUBJID = factor(rep(c("S1", "S2", "S3", "S4"), each = 3)),
    AVISIT = factor(
      rep(c("Week 4", "Week 8", "Week 12"), 4),
      levels = c("Week 4", "Week 8", "Week 12")
    ),
    TRT = factor(rep(c("Placebo", "Drug A"), each = 6)),
    CHG = c(1.1, 2.2, 3.3, 0.5, NA, NA, 1.0, 2.0, NA, 1.5, NA, 2.5),
    BASE = rep(c(10, 12, 11, 9), each = 3),
    STRATA = factor(rep(c("A", "B", "A", "B"), each = 3)),
    stringsAsFactors = FALSE
  )
}

make_test_vars <- function() {
  rbmi::set_vars(
    subjid = "USUBJID",
    visit = "AVISIT",
    group = "TRT",
    outcome = "CHG",
    covariates = c("BASE", "STRATA")
  )
}
```

**Teardown Pattern:**
- No explicit teardown blocks needed (testthat isolation)
- Fixture functions ensure fresh data for each test

**Assertion Pattern:**
- Descriptive test names explain what's being tested
- Multiple assertions per test block
- Clear expectation statements with context

## Mocking

**Framework:** None detected; tests use real objects

**Patterns:**
- Real data fixtures created via `make_test_data()`, `make_integration_data()`
- Real rbmi objects used in integration tests
- Package dependencies mocked via `testthat::skip_if_not_installed()`

**Example from `test-analysis_utils.R`:**
```r
testthat::test_that("gcomp_responder runs, drops visit from model terms, returns stable shape", {
  testthat::skip_if_not_installed("beeca")
  set.seed(1)
  dat <- data.frame(
    Y = rbinom(160, 1, 0.45),
    TRT = factor(
      sample(c("Placebo", "Drug"), 160, TRUE),
      levels = c("Placebo", "Drug")
    ),
    BASE = rnorm(160),
    VIS = sample(c("W4", "W8"), 160, TRUE)
  )
  vars <- list(
    outcome = "Y",
    group = "TRT",
    covariates = c("BASE", "TRT:BASE", "VIS"),
    visit = "VIS"
  )

  out <- gcomp_responder(
    data = dat,
    vars = vars,
    reference_levels = "Placebo",
    var_method = "Ge",
    type = "HC0",
    contrast = "diff"
  )

  # Assertions verify real behavior
  expect_true(any(grepl("^trt_", names(out))))
  expect_true(any(grepl("^lsm_", names(out))))
})
```

**What to Mock:**
- Not typical; prefer real objects and skip tests if dependencies unavailable

**What NOT to Mock:**
- rbmi objects (real imputation objects used)
- data.frame structures (fixtures created fresh)
- Dependencies optional via skip_if_not_installed()

## Fixtures and Factories

**Test Data:**
```r
# Factories for reusable test data (from test-integration.R)
make_integration_data <- function() {
  set.seed(42)
  n_subjects <- 50
  visits <- c("Week 4", "Week 8", "Week 12")

  data <- expand.grid(
    USUBJID = sprintf("S%03d", 1:n_subjects),
    AVISIT = visits,
    stringsAsFactors = FALSE
  )

  # ... populate with realistic patterns ...
  # ~20% dropout (monotone missing)
  # Treatment effect in outcome
  # Binary responder flag

  data
}
```

**Location:**
- Within each test file at the top
- Helper functions: `make_test_data()`, `make_test_vars()`, `make_integration_data()`, `make_integration_vars()`
- Example file: `tests/testthat/test-data_helpers.R` (lines 1-26)

## Coverage

**Requirements:** Not enforced (no codecov config detected)

**View Coverage:**
```bash
# From R console
devtools::test_coverage()

# Or via GitHub Actions if configured (not visible in current state)
```

**Current Coverage:**
- 159 test blocks across all modules
- Comprehensive coverage of:
  - Input validation (error cases)
  - Normal operation (happy paths)
  - Edge cases (boundary values, NA handling)
  - Integration workflows (multi-function pipelines)

## Test Types

**Unit Tests:**
- Scope: Individual functions
- Approach: Test each function with multiple scenarios
- Example: `test-formatting.R` tests `format_pvalue()`, `format_estimate()`, `format_results_table()` in isolation
- Coverage: Input validation, output format, edge cases (NA, boundary values)

```r
test_that("format_pvalue handles NA values", {
  expect_true(is.na(format_pvalue(NA_real_)))

  pvals <- c(0.05, NA, 0.001)
  result <- format_pvalue(pvals)
  expect_true(is.na(result[2]))
  expect_equal(result[1], "0.050")
})
```

**Integration Tests:**
- Scope: Multi-function workflows and complete analysis pipelines
- File: `tests/testthat/test-integration.R` (419 lines, 15+ test blocks)
- Approach: Test realistic clinical trial analysis workflows
- Coverage: Data validation → analysis → pooling → tidying pipelines

```r
test_that("Complete ANCOVA workflow: validate -> analyse -> pool -> tidy", {
  data("ADMI")
  set.seed(123)

  # Step 1: Prepare data
  ADMI$TRT <- factor(ADMI$TRT, levels = c("Placebo", "Drug A"))
  ADMI$USUBJID <- factor(ADMI$USUBJID)
  ADMI$AVISIT <- factor(ADMI$AVISIT)

  vars <- rbmi::set_vars(...)
  method <- rbmi::method_bayes(...)

  # Step 2: Validate
  expect_true(validate_data(ADMI, vars))

  # Step 3: Analyse
  ana_obj <- analyse_mi_data(data = ADMI, vars = vars, method = method, fun = rbmi::ancova)

  # Step 4: Pool
  pool_obj <- rbmi::pool(ana_obj)

  # Step 5: Tidy
  result <- tidy_pool_obj(pool_obj)

  # Assertions on final output
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("parameter", "est", "lci", "uci", "pval") %in% names(result)))
})
```

**E2E Tests:**
- Not explicitly labeled; integration tests serve this purpose
- Closest: `test-utils.R` line 1 test uses full rbmi workflow with real imputation

## Common Patterns

**Async Testing:**
Not applicable (R is single-threaded in tests)

**Error Testing:**
```r
# From test-formatting.R
test_that("format_pvalue validates inputs", {
  expect_error(format_pvalue("0.05"), "numeric")
  expect_error(format_pvalue(0.05, digits = -1), "positive")
  expect_error(format_pvalue(0.05, threshold = 0), "positive")
})

# From test-data_helpers.R
test_that("validate_data errors when data is not a data.frame", {
  vars <- make_test_vars()
  expect_error(validate_data(list(a = 1), vars), "must be a data.frame")
})

test_that("validate_data collects multiple issues in one error", {
  dat <- make_test_data()
  dat$CHG <- as.character(dat$CHG)
  dat$BASE[1] <- NA
  vars <- make_test_vars()
  err <- tryCatch(validate_data(dat, vars), error = function(e) e$message)
  expect_true(grepl("CHG", err))
  expect_true(grepl("BASE", err))
})
```

**Vector/Length Testing:**
```r
test_that("format_estimate handles vector input", {
  result <- format_estimate(
    estimate = c(-2.5, -1.8),
    lower = c(-4.0, -3.2),
    upper = c(-1.0, -0.4)
  )

  expect_length(result, 2)
  expect_equal(result[1], "-2.50 (-4.00, -1.00)")
  expect_equal(result[2], "-1.80 (-3.20, -0.40)")
})
```

**Conditional Testing:**
```r
# From test-analysis_utils.R
test_that("as_simple_formula2 builds intended formula", {
  # ...
  intercept_ok <- FALSE
  out <- try(as_simple_formula2("Y", character(0)), silent = TRUE)
  if (!inherits(out, "try-error")) {
    expect_identical(gsub("\\s+", "", deparse(out)), "Y~1")
    intercept_ok <- TRUE
  }
  if (!intercept_ok) {
    expect_error(as_simple_formula2("Y", character(0)))
  }
})
```

## Test Assertion Usage

**Most Common Assertions:**
- `expect_equal()`: 102 uses (equality of values)
- `expect_true()`: 94 uses (boolean conditions)
- `expect_error()`: 64 uses (error conditions)
- `expect_s3_class()`: 26 uses (class checking)
- `expect_type()`: 9 uses (type checking)

**Less Common:**
- `expect_false()`: 7 uses
- `expect_identical()`: 6 uses
- `expect_warning()`: 4 uses
- `expect_named()`: 4 uses
- `expect_match()`: 3 uses (regex matching)

## Test Data Examples

From `test-integration.R`, realistic clinical trial data factory:

```r
make_integration_data <- function() {
  set.seed(42)
  n_subjects <- 50
  visits <- c("Week 4", "Week 8", "Week 12")

  # Create all subject-visit combinations
  data <- expand.grid(
    USUBJID = sprintf("S%03d", 1:n_subjects),
    AVISIT = visits,
    stringsAsFactors = FALSE
  )

  # Factor conversion
  data$USUBJID <- factor(data$USUBJID)
  data$AVISIT <- factor(data$AVISIT, levels = visits)
  data$TRT <- factor(
    rep(c("Placebo", "Drug A"), each = n_subjects / 2)[as.integer(factor(data$USUBJID))],
    levels = c("Placebo", "Drug A")
  )

  # Outcome with structure (treatment effect + visit effect + noise)
  data$BASE <- rep(rnorm(n_subjects, 50, 10), each = length(visits))
  trt_effect <- ifelse(data$TRT == "Drug A", -5, 0)
  visit_effect <- as.integer(data$AVISIT) * 2
  data$CHG <- data$BASE * 0.1 + trt_effect + visit_effect + rnorm(nrow(data), 0, 5)

  # Monotone missing (~20% dropout)
  dropout_subjects <- sample(unique(data$USUBJID), n_subjects * 0.2)
  dropout_visits <- sample(2:3, length(dropout_subjects), replace = TRUE)
  for (i in seq_along(dropout_subjects)) {
    subj <- dropout_subjects[i]
    drop_visit <- dropout_visits[i]
    mask <- data$USUBJID == subj & as.integer(data$AVISIT) >= drop_visit
    data$CHG[mask] <- NA
  }

  # Binary responder
  data$RESP <- as.integer(data$CHG < 0)
  data$RESP[is.na(data$CHG)] <- NA

  data
}
```

---

*Testing analysis: 2026-02-07*
