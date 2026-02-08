# =============================================================================
# Tests for efficacy_table() function
# =============================================================================

# Helper to build a proper mock pool object matching rbmi internal structure
# (per decision 01-01-D3, copied from test-pool_methods.R for self-containment)
make_mock_pool <- function() {
  pool_obj <- list(
    pars = list(
      trt_Week4 = list(est = -2.5, se = 0.8, ci = c(-4.1, -0.9), pvalue = 0.002),
      lsm_ref_Week4 = list(est = 10.0, se = 0.5, ci = c(9.0, 11.0), pvalue = NA),
      lsm_alt_Week4 = list(est = 7.5, se = 0.6, ci = c(6.3, 8.7), pvalue = NA),
      trt_Week8 = list(est = -1.0, se = 1.2, ci = c(-3.4, 1.4), pvalue = 0.42)
    ),
    conf.level = 0.95,
    alternative = "two.sided",
    N = 100,
    method = "rubin"
  )
  class(pool_obj) <- "pool"
  pool_obj
}


# --- Dependency guard ---

test_that("efficacy_table errors when gt is not available", {
  skip_if_not_installed("gt")

  mock_pool <- make_mock_pool()

  local_mocked_bindings(is_gt_available = function() FALSE)

  expect_error(
    efficacy_table(mock_pool),
    class = "rbmiUtils_error_dependency"
  )
})


# --- Input validation ---

test_that("efficacy_table errors on non-pool input", {
  skip_if_not_installed("gt")

  expect_error(
    efficacy_table("not_a_pool"),
    class = "rbmiUtils_error_validation"
  )

  expect_error(
    efficacy_table(42),
    class = "rbmiUtils_error_validation"
  )
})


# --- Returns gt_tbl object ---

test_that("efficacy_table returns gt_tbl object", {
  skip_if_not_installed("gt")

  mock_pool <- make_mock_pool()
  tbl <- efficacy_table(mock_pool)
  expect_s3_class(tbl, "gt_tbl")
})


# --- Visit labels ---

test_that("table contains cleaned visit labels", {
  skip_if_not_installed("gt")

  mock_pool <- make_mock_pool()
  tbl <- efficacy_table(mock_pool)
  html <- gt::as_raw_html(tbl)

  # Visit names should be cleaned: "Week4" -> "Week 4"
  expect_true(grepl("Week 4", html, fixed = TRUE))
  expect_true(grepl("Week 8", html, fixed = TRUE))
})


# --- Row labels ---

test_that("table contains correct row labels", {
  skip_if_not_installed("gt")

  html <- gt::as_raw_html(efficacy_table(make_mock_pool()))

  expect_true(grepl("LS Mean", html))
  expect_true(grepl("Treatment Difference", html))
})


# --- Custom arm labels ---

test_that("custom arm labels appear in table", {
  skip_if_not_installed("gt")

  tbl <- efficacy_table(
    make_mock_pool(),
    arm_labels = c(ref = "Placebo", alt = "Drug A")
  )
  html <- gt::as_raw_html(tbl)

  expect_true(grepl("Placebo", html))
  expect_true(grepl("Drug A", html))
})


# --- Default arm labels ---

test_that("default arm labels when not provided", {
  skip_if_not_installed("gt")

  html <- gt::as_raw_html(efficacy_table(make_mock_pool()))

  expect_true(grepl("Reference", html))
  expect_true(grepl("Treatment", html))
})


# --- Footnotes ---

test_that("footnotes are present", {
  skip_if_not_installed("gt")

  html <- gt::as_raw_html(efficacy_table(make_mock_pool()))

  # Pooling method
  expect_true(grepl("rubin", html, ignore.case = TRUE))
  # Number of imputations
  expect_true(grepl("100", html))
  # Confidence level
  expect_true(grepl("95%", html))
})


# --- Title and subtitle ---

test_that("title and subtitle appear", {
  skip_if_not_installed("gt")

  tbl <- efficacy_table(
    make_mock_pool(),
    title = "Table 14.2.1",
    subtitle = "ANCOVA Model"
  )
  html <- gt::as_raw_html(tbl)

  expect_true(grepl("Table 14.2.1", html))
  expect_true(grepl("ANCOVA Model", html))
})


# --- P-value formatting ---

test_that("p-value formatting: em dash for LSM, formatted for trt", {
  skip_if_not_installed("gt")

  html <- gt::as_raw_html(efficacy_table(make_mock_pool()))

  # Em dash present (for LSM rows with NA p-values)
  # Could be encoded as UTF-8, HTML entity, or raw character
  has_em_dash <- grepl("\u2014", html) ||
    grepl("&mdash;", html) ||
    grepl("&#8212;", html) ||
    grepl("\xe2\x80\x94", html, useBytes = TRUE)
  expect_true(has_em_dash)

  # Formatted p-value present (0.002 from trt_Week4, 0.420 from trt_Week8)
  has_pval <- grepl("0.002", html) || grepl("0.420", html)
  expect_true(has_pval)
})


# --- Digits parameter ---

test_that("digits parameter controls formatting", {
  skip_if_not_installed("gt")

  tbl3 <- efficacy_table(make_mock_pool(), digits = 3)
  html3 <- gt::as_raw_html(tbl3)

  # With 3 decimal places, the CI should show 3 decimal formatting
  # e.g., -4.100 or 9.000
  expect_true(grepl("\\d+\\.\\d{3}", html3))
})


# --- Custom CI level ---

test_that("custom ci_level label appears in table", {
  skip_if_not_installed("gt")

  mock_pool <- make_mock_pool()
  mock_pool$conf.level <- 0.90

  html <- gt::as_raw_html(efficacy_table(mock_pool))

  # 90% CI should be the column label and footnote
  expect_true(grepl("90%", html))
})


# --- CI level from argument overrides pool object ---

test_that("ci_level argument overrides pool_obj$conf.level", {
  skip_if_not_installed("gt")

  mock_pool <- make_mock_pool()
  # pool_obj has conf.level = 0.95
  html <- gt::as_raw_html(efficacy_table(mock_pool, ci_level = 0.99))

  expect_true(grepl("99%", html))
})
