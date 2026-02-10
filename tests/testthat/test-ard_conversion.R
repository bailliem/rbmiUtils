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


# =============================================================================
# Tests for compute_rubin_diagnostics() internal function
# =============================================================================

test_that("compute_rubin_diagnostics returns correct values for known inputs", {
  # Test case 1: Known values with hand-calculated expected results
  M <- 5L
  ests <- c(1.0, 1.2, 0.8, 1.1, 0.9)
  ses <- c(0.5, 0.6, 0.4, 0.55, 0.45)
  v_com <- 100

  result <- compute_rubin_diagnostics(ests, ses, v_com, M)

  # Should return a named list

  expect_type(result, "list")
  expect_named(result, c("var_w", "var_b", "var_t", "lambda", "riv",
                          "df_adj", "dfcom", "fmi", "re"))

  # Verify intermediates
  expect_equal(result$var_w, 0.255, tolerance = 1e-10)
  expect_equal(result$var_b, 0.025, tolerance = 1e-10)
  expect_equal(result$var_t, 0.285, tolerance = 1e-10)

  # Verify derived quantities
  expect_equal(result$lambda, 0.105263157894737, tolerance = 1e-10)
  expect_equal(result$riv, 0.117647058823529, tolerance = 1e-10)
  expect_equal(result$df_adj, 70.582240254527292, tolerance = 1e-6)
  expect_equal(result$dfcom, 100)
  expect_equal(result$fmi, 0.129582527324379, tolerance = 1e-6)
  expect_equal(result$re, 0.974738192312120, tolerance = 1e-6)
})


test_that("compute_rubin_diagnostics handles lambda=0 (identical estimates)", {
  # Test case 2: All estimates identical -> var_b = 0, lambda = 0
  M <- 5L
  ests <- rep(1.0, 5)
  ses <- c(0.5, 0.6, 0.4, 0.55, 0.45)
  v_com <- 100

  result <- compute_rubin_diagnostics(ests, ses, v_com, M)

  expect_equal(result$var_b, 0)
  expect_equal(result$lambda, 0)
  expect_equal(result$riv, 0)

  # df_adj should be v_obs when lambda=0
  expected_v_obs <- ((v_com + 1) / (v_com + 3)) * v_com * (1 - 0)
  expect_equal(result$df_adj, expected_v_obs, tolerance = 1e-10)

  # fmi should be near 0 (small positive due to df adjustment term)
  expect_true(result$fmi >= 0)
  expect_true(result$fmi < 0.1)

  # re should be near 1
  expect_true(result$re > 0.99)
  expect_true(result$re <= 1.0)
})


test_that("compute_rubin_diagnostics handles v_com=Inf (large-sample)", {
  # Test case 3: v_com = Inf with lambda > 0
  M <- 5L
  ests <- c(1.0, 1.2, 0.8, 1.1, 0.9)
  ses <- c(0.5, 0.6, 0.4, 0.55, 0.45)
  v_com <- Inf

  result <- compute_rubin_diagnostics(ests, ses, v_com, M)

  # df_adj should be v_old = (M-1)/lambda^2 when v_com is Inf and lambda > 0
  expected_v_old <- (M - 1) / result$lambda^2
  expect_equal(result$df_adj, expected_v_old, tolerance = 1e-6)

  # dfcom pass-through
  expect_equal(result$dfcom, Inf)

  # fmi and re should be valid numeric
  expect_true(is.finite(result$fmi))
  expect_true(is.finite(result$re))
  expect_equal(result$fmi, 0.110179294389821, tolerance = 1e-6)
  expect_equal(result$re, 0.978439250749816, tolerance = 1e-6)
})


test_that("compute_rubin_diagnostics handles v_com=Inf and lambda=0", {
  # Edge case: v_com = Inf AND var_b = 0 -> df_adj = Inf
  M <- 5L
  ests <- rep(1.0, 5)
  ses <- c(0.5, 0.6, 0.4, 0.55, 0.45)
  v_com <- Inf

  result <- compute_rubin_diagnostics(ests, ses, v_com, M)

  expect_equal(result$var_b, 0)
  expect_equal(result$lambda, 0)
  expect_true(is.infinite(result$df_adj))
  expect_equal(result$dfcom, Inf)
})


test_that("compute_rubin_diagnostics handles all NA SEs", {
  # Test case 4: All SEs are NA -> all diagnostics should be NA
  M <- 5L
  ests <- c(1.0, 1.2, 0.8, 1.1, 0.9)
  ses <- rep(NA_real_, 5)
  v_com <- 100

  result <- compute_rubin_diagnostics(ests, ses, v_com, M)

  # All values should be NA
  expect_true(is.na(result$var_w))
  expect_true(is.na(result$lambda))
  expect_true(is.na(result$riv))
  expect_true(is.na(result$df_adj))
  expect_true(is.na(result$fmi))
  expect_true(is.na(result$re))

  # dfcom should still be passed through
  expect_equal(result$dfcom, v_com)
})


test_that("compute_rubin_diagnostics: fmi is distinct from lambda", {
  # Test case 5: Verify fmi != lambda (adjusted FMI differs from lambda)
  M <- 5L
  ests <- c(1.0, 1.2, 0.8, 1.1, 0.9)
  ses <- c(0.5, 0.6, 0.4, 0.55, 0.45)
  v_com <- 100

  result <- compute_rubin_diagnostics(ests, ses, v_com, M)

  # fmi and lambda must both be numeric
  expect_true(is.numeric(result$fmi))
  expect_true(is.numeric(result$lambda))

  # They must be different (adjusted FMI accounts for finite df)
  expect_false(isTRUE(all.equal(result$fmi, result$lambda)),
               info = "FMI and lambda should be distinct quantities")

  # fmi should be greater than lambda (adjusted FMI >= lambda in general)
  expect_true(result$fmi > result$lambda)
})


test_that("compute_rubin_diagnostics handles NA v_com", {
  # Edge case: v_com is NA
  M <- 5L
  ests <- c(1.0, 1.2, 0.8, 1.1, 0.9)
  ses <- c(0.5, 0.6, 0.4, 0.55, 0.45)
  v_com <- NA_real_

  result <- compute_rubin_diagnostics(ests, ses, v_com, M)

  # var_w, var_b, var_t, lambda, riv should still be computable
  expect_true(is.finite(result$var_w))
  expect_true(is.finite(result$var_b))
  expect_true(is.finite(result$var_t))
  expect_true(is.finite(result$lambda))
  expect_true(is.finite(result$riv))

  # df_adj, fmi, re should be NA when v_com is NA
  expect_true(is.na(result$df_adj))
  expect_true(is.na(result$fmi))
  expect_true(is.na(result$re))

  # dfcom should be NA (pass-through)
  expect_true(is.na(result$dfcom))
})
