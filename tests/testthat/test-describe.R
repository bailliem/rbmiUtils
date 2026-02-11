# tests/testthat/test-describe.R
# Tests for describe_draws() and print.describe_draws()

# ---------------------------------------------------------------------------
# Mock object helpers
# ---------------------------------------------------------------------------

#' Create a minimal mock sample_single object
make_mock_sample <- function(failed = FALSE, n_ids = 5) {
  ids <- paste0("SUBJ", seq_len(n_ids))
  if (failed) {
    structure(
      list(
        failed = TRUE,
        ids = ids,
        ids_samp = ids,
        beta = NA,
        sigma = NA,
        theta = NA
      ),
      class = "sample_single"
    )
  } else {
    sigma_list <- list(grp1 = matrix(1, 2, 2))
    structure(
      list(
        failed = FALSE,
        ids = ids,
        ids_samp = ids,
        beta = c(intercept = 1.0, trt = -0.5),
        sigma = sigma_list,
        theta = c(1.0, 0.5, 0.2)
      ),
      class = "sample_single"
    )
  }
}

#' Create a mock draws object for condmean method
#' @param type "jackknife" or "bootstrap"
#' @param n_resampled Number of resampled draws (total samples = 1 + n_resampled)
#' @param n_failures Number of failed samples
make_mock_draws_condmean <- function(type = "jackknife",
                                     n_resampled = 10,
                                     n_failures = 0) {
  method <- structure(
    list(
      covariance = "us",
      threshold = 0.01,
      same_cov = TRUE,
      REML = TRUE,
      type = type
    ),
    class = c("method", "condmean")
  )

  n_total <- 1 + n_resampled
  samples <- lapply(seq_len(n_total), function(i) make_mock_sample())
  class(samples) <- "sample_list"

  structure(
    list(
      data = NULL,
      method = method,
      samples = samples,
      n_failures = n_failures,
      fit = NULL,
      formula = CHG ~ 1 + TRT + AVISIT + BASE
    ),
    class = c("condmean", "draws", "list")
  )
}

#' Create a mock draws object for bayes method (no stanfit)
#' @param n_samples Number of samples
#' @param n_failures Number of failed samples
make_mock_draws_bayes <- function(n_samples = 100, n_failures = 0) {
  method <- structure(
    list(
      covariance = "us",
      same_cov = FALSE,
      n_samples = n_samples,
      prior_cov = NULL,
      control = list(
        warmup = 200,
        thin = 1,
        chains = 4,
        seed = 123
      )
    ),
    class = c("method", "bayes")
  )

  samples <- lapply(seq_len(n_samples), function(i) make_mock_sample())
  class(samples) <- "sample_list"

  structure(
    list(
      data = NULL,
      method = method,
      samples = samples,
      n_failures = n_failures,
      fit = NULL,
      formula = CHG ~ 1 + TRT + AVISIT + BASE
    ),
    class = c("random", "draws", "list")
  )
}

#' Create a mock draws object for approxbayes method
#' @param n_samples Number of samples
#' @param n_failures Number of failed samples
make_mock_draws_approxbayes <- function(n_samples = 100, n_failures = 0) {
  method <- structure(
    list(
      covariance = "us",
      threshold = 0.01,
      same_cov = TRUE,
      REML = TRUE,
      n_samples = n_samples
    ),
    class = c("method", "approxbayes")
  )

  samples <- lapply(seq_len(n_samples), function(i) make_mock_sample())
  class(samples) <- "sample_list"

  structure(
    list(
      data = NULL,
      method = method,
      samples = samples,
      n_failures = n_failures,
      fit = NULL,
      formula = CHG ~ 1 + TRT + AVISIT + BASE
    ),
    class = c("random", "draws", "list")
  )
}


# ===========================================================================
# Test: Input validation
# ===========================================================================

test_that("describe_draws rejects non-draws objects", {
  expect_error(
    describe_draws(list(a = 1)),
    class = "rbmiUtils_error_type"
  )

  expect_error(
    describe_draws("not a draws"),
    class = "rbmiUtils_error_type"
  )

  expect_error(
    describe_draws(NULL),
    class = "rbmiUtils_error_type"
  )
})


# ===========================================================================
# Test: S3 class of result
# ===========================================================================

test_that("describe_draws returns S3 class c('describe_draws', 'list')", {
  draws_obj <- make_mock_draws_condmean()
  desc <- describe_draws(draws_obj)

  expect_s3_class(desc, "describe_draws")
  expect_true(inherits(desc, "list"))
  expect_equal(class(desc), c("describe_draws", "list"))
})


# ===========================================================================
# Test: Condmean jackknife
# ===========================================================================

test_that("describe_draws extracts correct info for condmean jackknife", {
  draws_obj <- make_mock_draws_condmean(type = "jackknife", n_resampled = 10)
  desc <- describe_draws(draws_obj)

  expect_equal(desc$method, "Conditional Mean (jackknife)")
  expect_equal(desc$method_class, "condmean")
  expect_equal(desc$n_samples, 11)
  expect_equal(desc$n_failures, 0)
  expect_equal(desc$covariance, "us")
  expect_true(desc$same_cov)
  expect_match(desc$formula, "CHG")

  # Condmean-specific fields

  expect_equal(desc$condmean_type, "jackknife")
  expect_equal(desc$n_primary, 1)
  expect_equal(desc$n_resampled, 10)
})


# ===========================================================================
# Test: Condmean bootstrap
# ===========================================================================

test_that("describe_draws returns correct condmean_type for bootstrap", {
  draws_obj <- make_mock_draws_condmean(type = "bootstrap", n_resampled = 5)
  desc <- describe_draws(draws_obj)

  expect_equal(desc$method, "Conditional Mean (bootstrap)")
  expect_equal(desc$condmean_type, "bootstrap")
  expect_equal(desc$n_primary, 1)
  expect_equal(desc$n_resampled, 5)
  expect_equal(desc$n_samples, 6)
})


# ===========================================================================
# Test: Condmean with 0 resampled
# ===========================================================================

test_that("describe_draws handles condmean with 0 resampled samples", {
  draws_obj <- make_mock_draws_condmean(type = "jackknife", n_resampled = 0)
  desc <- describe_draws(draws_obj)

  expect_equal(desc$n_samples, 1)
  expect_equal(desc$n_primary, 1)

  expect_equal(desc$n_resampled, 0)
})


# ===========================================================================
# Test: Bayesian without rstan / without stanfit
# ===========================================================================

test_that("describe_draws returns NULL mcmc when fit is NULL (bayes)", {
  draws_obj <- make_mock_draws_bayes(n_samples = 50)
  desc <- describe_draws(draws_obj)

  expect_equal(desc$method, "Bayesian (MCMC via Stan)")
  expect_equal(desc$method_class, "bayes")
  expect_equal(desc$n_samples, 50)
  expect_null(desc$mcmc)

  # Bayes-specific fields
  expect_equal(desc$bayes_control$warmup, 200)
  expect_equal(desc$bayes_control$thin, 1)
  expect_equal(desc$bayes_control$chains, 4)
  expect_equal(desc$bayes_control$seed, 123)
})


# ===========================================================================
# Test: Bayesian with stanfit (requires rstan)
# ===========================================================================

test_that("describe_draws extracts MCMC diagnostics from stanfit", {
  skip_if_not_installed("rstan")

  # Create a mock stanfit-like object that rstan::summary can handle
  # This test uses a real stanfit if available from a minimal model
  # For now, we create a mock that simulates the describe_draws mcmc extraction


  # Build a mock draws with a fake stanfit to test the extraction logic
  # We test the actual stanfit integration separately
  draws_obj <- make_mock_draws_bayes(n_samples = 50)

  # Simulate a stanfit-like object: we need draws_obj$fit to be "stanfit" class
  # and rstan::summary to work on it. Instead, we test with mocking the path.
  # If rstan is installed, try building a real tiny model
  skip("Stanfit mock requires real Stan compilation - tested via integration")
})


# ===========================================================================
# Test: approxbayes (no MCMC diagnostics)
# ===========================================================================

test_that("describe_draws returns NULL mcmc for approxbayes", {
  draws_obj <- make_mock_draws_approxbayes(n_samples = 80)
  desc <- describe_draws(draws_obj)

  expect_equal(desc$method, "Approximate Bayesian")
  expect_equal(desc$method_class, "approxbayes")
  expect_equal(desc$n_samples, 80)
  expect_null(desc$mcmc)

  # approxbayes should not have condmean-specific fields
  expect_null(desc$condmean_type)
  expect_null(desc$n_primary)
  expect_null(desc$n_resampled)

  # approxbayes should not have bayes_control
  expect_null(desc$bayes_control)
})


# ===========================================================================
# Test: Failures count
# ===========================================================================

test_that("describe_draws reports failure count correctly", {
  draws_obj <- make_mock_draws_condmean(n_resampled = 10, n_failures = 3)
  desc <- describe_draws(draws_obj)

  expect_equal(desc$n_failures, 3)
})


# ===========================================================================
# Test: print method returns invisible(x)
# ===========================================================================

test_that("print.describe_draws returns invisible(x) and produces output", {
  draws_obj <- make_mock_draws_condmean(type = "jackknife", n_resampled = 10)
  desc <- describe_draws(draws_obj)

  # Capture cli output (writes to stderr/message connection)
  out <- capture.output(result <- print(desc), type = "message")

  # Returns invisible(x) - same object
  expect_identical(result, desc)

  # Produces some output
  expect_true(length(out) > 0)
})


# ===========================================================================
# Test: print method shows correct content for condmean
# ===========================================================================

test_that("print.describe_draws shows 1 + N format for condmean", {
  draws_obj <- make_mock_draws_condmean(type = "jackknife", n_resampled = 10)
  desc <- describe_draws(draws_obj)

  out <- paste(capture.output(print(desc), type = "message"), collapse = "\n")

  # Should contain "1 + 10" format
  expect_match(out, "1 \\+ 10", fixed = FALSE)
  expect_match(out, "jackknife", ignore.case = TRUE)
})


# ===========================================================================
# Test: print method shows Bayesian info
# ===========================================================================

test_that("print.describe_draws shows Bayesian method info", {
  draws_obj <- make_mock_draws_bayes(n_samples = 50)
  desc <- describe_draws(draws_obj)

  out <- paste(capture.output(print(desc), type = "message"), collapse = "\n")

  expect_match(out, "Bayesian", ignore.case = TRUE)
  expect_match(out, "50")
})


# ===========================================================================
# Test: Formula extraction
# ===========================================================================

test_that("describe_draws extracts formula as character string", {
  draws_obj <- make_mock_draws_condmean()
  desc <- describe_draws(draws_obj)

  expect_type(desc$formula, "character")
  expect_match(desc$formula, "CHG")
  expect_match(desc$formula, "TRT")
  expect_match(desc$formula, "BASE")
})


# ===========================================================================
# Test: same_cov field
# ===========================================================================

test_that("describe_draws captures same_cov correctly", {
  # condmean has same_cov = TRUE by default in mock
  draws_obj <- make_mock_draws_condmean()
  desc <- describe_draws(draws_obj)
  expect_true(desc$same_cov)

  # bayes mock has same_cov = FALSE
  draws_obj2 <- make_mock_draws_bayes()
  desc2 <- describe_draws(draws_obj2)
  expect_false(desc2$same_cov)
})
