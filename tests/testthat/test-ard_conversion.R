# =============================================================================
# Tests for pool_to_ard() ARD conversion
# =============================================================================

# Helper to build a proper mock pool object matching rbmi internal structure
# (per decision 01-01-D3)
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


test_that("pool_to_ard returns valid ARD structure", {
  skip_if_not_installed("cards")

  mock_pool <- make_mock_pool()
  ard <- pool_to_ard(mock_pool)

  # Result has class "card"
  expect_s3_class(ard, "card")

  # Result has required columns
  required_cols <- c(
    "group1", "group1_level",
    "group2", "group2_level",
    "group3", "group3_level",
    "variable", "variable_level",
    "stat_name", "stat_label",
    "stat", "warning", "error"
  )
  for (col in required_cols) {
    expect_true(col %in% names(ard), info = paste("Missing column:", col))
  }

  # cards::check_ard_structure() does not error
  expect_no_error(cards::check_ard_structure(ard))
})


test_that("pool_to_ard includes all five statistics per parameter", {
  skip_if_not_installed("cards")

  mock_pool <- make_mock_pool()
  ard <- pool_to_ard(mock_pool)

  expected_stats <- c("estimate", "std.error", "conf.low", "conf.high", "p.value")
  unique_vars <- unique(ard$variable)

  for (v in unique_vars) {
    var_ard <- ard[ard$variable == v, ]
    var_stat_names <- var_ard$stat_name

    # Each variable should have the 5 statistics (plus method)
    for (s in expected_stats) {
      expect_true(
        s %in% var_stat_names,
        info = paste("Variable", v, "missing stat_name:", s)
      )
    }
    # Also has method row
    expect_true("method" %in% var_stat_names,
                info = paste("Variable", v, "missing method row"))
  }
})


test_that("pool_to_ard preserves grouping columns (ARD-03)", {
  skip_if_not_installed("cards")

  mock_pool <- make_mock_pool()
  ard <- pool_to_ard(mock_pool)

  # group1 is always "visit"
  expect_true(all(ard$group1 == "visit"))

  # group2 is always "parameter_type"
  expect_true(all(ard$group2 == "parameter_type"))

  # group3 is always "lsm_type"
  expect_true(all(ard$group3 == "lsm_type"))

  # group1_level values should match visits from tidy_pool_obj
  tidy_df <- tidy_pool_obj(mock_pool)
  expected_visits <- sort(unique(tidy_df$visit))
  ard_visits <- sort(unique(unlist(ard$group1_level)))
  expect_equal(ard_visits, expected_visits)

  # parameter_type values should include both "trt" and "lsm"
  ard_ptypes <- unique(unlist(ard$group2_level))
  expect_true("trt" %in% ard_ptypes)
  expect_true("lsm" %in% ard_ptypes)
})


test_that("pool_to_ard stat column is list column", {
  skip_if_not_installed("cards")

  mock_pool <- make_mock_pool()
  ard <- pool_to_ard(mock_pool)

  # stat must be a list column (Pitfall 2)
  expect_true(is.list(ard$stat))

  # group level columns must be list columns (Pitfall 3)
  expect_true(is.list(ard$group1_level))
  expect_true(is.list(ard$group2_level))
  expect_true(is.list(ard$group3_level))
  expect_true(is.list(ard$variable_level))

  # warning and error must be list columns
  expect_true(is.list(ard$warning))
  expect_true(is.list(ard$error))

  # fmt_fun must be a list column
  expect_true(is.list(ard$fmt_fun))
})


test_that("pool_to_ard handles NA p-values", {
  skip_if_not_installed("cards")

  mock_pool <- make_mock_pool()
  ard <- pool_to_ard(mock_pool)

  # lsm_ref_Week4 has pvalue = NA
  lsm_ref_rows <- ard[ard$variable == "lsm_ref_Week4", ]
  pval_row <- lsm_ref_rows[lsm_ref_rows$stat_name == "p.value", ]

  # stat should contain NA (not NULL, not error)
  expect_true(is.na(pval_row$stat[[1]]))

  # lsm_alt_Week4 also has pvalue = NA
  lsm_alt_rows <- ard[ard$variable == "lsm_alt_Week4", ]
  pval_row2 <- lsm_alt_rows[lsm_alt_rows$stat_name == "p.value", ]
  expect_true(is.na(pval_row2$stat[[1]]))
})


test_that("pool_to_ard errors without cards package", {
  # Mock the internal helper that checks for cards availability
  local_mocked_bindings(
    is_cards_available = function() FALSE
  )

  mock_pool <- make_mock_pool()
  expect_error(
    pool_to_ard(mock_pool),
    class = "rbmiUtils_error_dependency"
  )
})


test_that("pool_to_ard validates input class", {
  skip_if_not_installed("cards")

  # Plain list without pool class
  bad_input <- list(pars = list(), method = "rubin")
  expect_error(
    pool_to_ard(bad_input),
    class = "rbmiUtils_error_validation"
  )

  # Numeric input
  expect_error(
    pool_to_ard(42),
    class = "rbmiUtils_error_validation"
  )

  # NULL input
  expect_error(
    pool_to_ard(NULL),
    class = "rbmiUtils_error_validation"
  )
})


test_that("pool_to_ard numeric values match tidy_pool_obj", {
  skip_if_not_installed("cards")

  mock_pool <- make_mock_pool()
  ard <- pool_to_ard(mock_pool)
  tidy_df <- tidy_pool_obj(mock_pool)

  # For each parameter in tidy_df, check that ARD estimate matches

  for (i in seq_len(nrow(tidy_df))) {
    param <- tidy_df$parameter[i]
    expected_est <- tidy_df$est[i]
    expected_se <- tidy_df$se[i]
    expected_lci <- tidy_df$lci[i]
    expected_uci <- tidy_df$uci[i]

    param_ard <- ard[ard$variable == param, ]

    # Estimate
    est_row <- param_ard[param_ard$stat_name == "estimate", ]
    expect_equal(est_row$stat[[1]], expected_est,
                 info = paste("Estimate mismatch for", param))

    # Std. Error
    se_row <- param_ard[param_ard$stat_name == "std.error", ]
    expect_equal(se_row$stat[[1]], expected_se,
                 info = paste("SE mismatch for", param))

    # Conf. Low
    lci_row <- param_ard[param_ard$stat_name == "conf.low", ]
    expect_equal(lci_row$stat[[1]], expected_lci,
                 info = paste("LCI mismatch for", param))

    # Conf. High
    uci_row <- param_ard[param_ard$stat_name == "conf.high", ]
    expect_equal(uci_row$stat[[1]], expected_uci,
                 info = paste("UCI mismatch for", param))
  }
})
