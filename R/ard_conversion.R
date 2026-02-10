#' Convert Pool Object to ARD Format
#'
#' Converts an rbmi pool object to the pharmaverse Analysis Results Dataset
#' (ARD) standard using the \pkg{cards} package. The ARD format is a
#' long-format data frame where each row represents a single statistic for a
#' given parameter, with grouping columns for visit, parameter type, and
#' least-squares-mean type.
#'
#' @param pool_obj A pooled analysis object of class `"pool"`, typically
#'   obtained from [rbmi::pool()] after calling [analyse_mi_data()].
#' @param conf.level Confidence level used for CI labels (e.g., `0.95` produces
#'   "95% CI Lower"). If `NULL` (the default), the value is taken from
#'   `pool_obj$conf.level`. If that is also `NULL`, defaults to `0.95`.
#'
#' @return A data frame of class `"card"` (ARD format) with grouping columns
#'   for visit (`group1`), parameter_type (`group2`), and lsm_type (`group3`).
#'   Each parameter produces rows for five statistics: estimate, std.error,
#'   conf.low, conf.high, and p.value, plus a method row.
#'
#' @details
#' The function works by:
#' 1. Tidying the pool object via [tidy_pool_obj()]
#' 2. Reshaping each parameter into long-format ARD rows (one row per statistic)
#' 3. Adding grouping columns (visit, parameter_type, lsm_type)
#' 4. Applying [cards::as_card()] and [cards::tidy_ard_column_order()] for
#'    standard ARD structure
#'
#' The resulting ARD passes [cards::check_ard_structure()] validation and is
#' suitable for downstream use with \pkg{gtsummary}.
#'
#' @seealso
#' * [rbmi::pool()] for creating pool objects
#' * [tidy_pool_obj()] for the underlying data transformation
#' * [cards::as_card()] and [cards::check_ard_structure()] for ARD validation
#'
#' @examples
#' \donttest{
#' # Requires the cards package
#' if (requireNamespace("cards", quietly = TRUE)) {
#'   # After running an rbmi analysis pipeline:
#'   # pool_obj <- rbmi::pool(analysis_obj)
#'   # ard <- pool_to_ard(pool_obj)
#'   # cards::check_ard_structure(ard)
#' }
#' }
#'
#' @export
pool_to_ard <- function(pool_obj, conf.level = NULL) {

  # --- Dependency check ---
  if (!is_cards_available()) {
    cli::cli_abort(
      c(
        "Package {.pkg cards} is required for ARD conversion.",
        "i" = "Install with {.code install.packages(\"cards\")}."
      ),
      class = c("rbmiUtils_error_dependency", "rbmiUtils_error")
    )
  }

  # --- Input validation ---
  if (!inherits(pool_obj, "pool")) {
    cli::cli_abort(
      "Input {.arg pool_obj} must be of class {.cls pool}, not {.cls {class(pool_obj)}}.",
      class = c("rbmiUtils_error_validation", "rbmiUtils_error")
    )
  }

  # --- Extract conf.level ---
  if (is.null(conf.level)) {
    conf.level <- pool_obj$conf.level %||% 0.95
  }

  # --- Tidy the pool object ---
  tidy_df <- tidy_pool_obj(pool_obj)

  # --- Define statistic names and labels ---
  stat_names <- c("estimate", "std.error", "conf.low", "conf.high", "p.value")
  stat_labels <- c(
    "Estimate",
    "Std. Error",
    paste0(conf.level * 100, "% CI Lower"),
    paste0(conf.level * 100, "% CI Upper"),
    "p-value"
  )
  n_stats <- length(stat_names)

  # --- Extract method label ---
  method_label <- if (!is.null(pool_obj$method)) pool_obj$method else "unknown"

  # --- Build long-format ARD rows ---
  rows <- lapply(seq_len(nrow(tidy_df)), function(i) {
    r <- tidy_df[i, ]

    lsm_val <- if (is.na(r$lsm_type)) NA_character_ else r$lsm_type

    # Number of rows = n_stats + 1 (for method row)
    n_total <- n_stats + 1L

    data.frame(
      group1       = rep("visit", n_total),
      group1_level = I(as.list(rep(r$visit, n_total))),
      group2       = rep("parameter_type", n_total),
      group2_level = I(as.list(rep(r$parameter_type, n_total))),
      group3       = rep("lsm_type", n_total),
      group3_level = I(as.list(rep(lsm_val, n_total))),
      variable       = rep(r$parameter, n_total),
      variable_level = I(as.list(rep(NA, n_total))),
      context        = rep("rbmi_pool", n_total),
      stat_name  = c(stat_names, "method"),
      stat_label = c(stat_labels, "Method"),
      stat       = I(list(r$est, r$se, r$lci, r$uci, r$pval, method_label)),
      fmt_fun    = I(as.list(rep(1L, n_total))),
      warning    = I(as.list(rep(list(NULL), n_total))),
      error      = I(as.list(rep(list(NULL), n_total))),
      stringsAsFactors = FALSE
    )
  })

  # --- Assemble and convert to ARD ---
  ard_df <- do.call(rbind, rows)
  cards::tidy_ard_column_order(cards::as_card(ard_df))
}


#' Check if cards package is available
#'
#' Internal helper extracted for testability. Wraps
#' `requireNamespace("cards", quietly = TRUE)`.
#'
#' @return Logical scalar.
#' @keywords internal
#' @noRd
is_cards_available <- function() {
  requireNamespace("cards", quietly = TRUE)
}


#' Compute Rubin's Rules Diagnostic Statistics
#'
#' Pure computational function implementing Rubin's rules variance
#' decomposition for multiple imputation diagnostics. Computes within-
#' and between-imputation variance, fraction of missing information (FMI),
#' lambda, relative increase in variance (RIV), Barnard-Rubin adjusted
#' degrees of freedom, and relative efficiency.
#'
#' @param ests Numeric vector of per-imputation point estimates.
#' @param ses Numeric vector of per-imputation standard errors.
#' @param v_com Numeric scalar for complete-data degrees of freedom.
#' @param M Integer, number of imputations.
#'
#' @return Named list with elements: `var_w` (within-imputation variance),
#'   `var_b` (between-imputation variance), `var_t` (total variance),
#'   `lambda` (proportion of variance due to missingness),
#'   `riv` (relative increase in variance),
#'   `df_adj` (Barnard-Rubin adjusted degrees of freedom),
#'   `dfcom` (complete-data df, pass-through of `v_com`),
#'   `fmi` (adjusted fraction of missing information),
#'   `re` (relative efficiency).
#'
#' @details
#' Formulas follow Rubin (1987) and Barnard & Rubin (1999), verified
#' against `rbmi:::rubin_rules` and `rbmi:::rubin_df` source code.
#' The adjusted FMI follows the mice package convention:
#' `fmi = (riv + 2/(df_adj + 3)) / (1 + riv)`.
#'
#' Edge case handling:
#' - All SEs are NA: returns list with all values NA (except dfcom)
#' - `var_b == 0` (lambda == 0): `df_adj = v_obs` (skip `v_old`)
#' - `v_com` is Inf and `var_b == 0`: `df_adj = Inf`
#' - `v_com` is Inf and `lambda > 0`: `df_adj = v_old = (M-1)/lambda^2`
#' - `v_com` is NA: `df_adj = NA`, `fmi = NA`, `re = NA`
#'
#' @keywords internal
#' @noRd
compute_rubin_diagnostics <- function(ests, ses, v_com, M) {

  # --- Handle all-NA SEs ---
  if (all(is.na(ses))) {
    return(list(
      var_w  = NA_real_,
      var_b  = var(ests),
      var_t  = NA_real_,
      lambda = NA_real_,
      riv    = NA_real_,
      df_adj = NA_real_,
      dfcom  = v_com,
      fmi    = NA_real_,
      re     = NA_real_
    ))
  }

  # --- Core Rubin's rules decomposition ---
  var_w <- mean(ses^2)                    # Within-imputation variance
  var_b <- var(ests)                      # Between-imputation variance
  var_t <- var_w + var_b + var_b / M      # Total variance

  # Lambda: proportion of total variance due to missingness
  lambda <- (1 + 1 / M) * var_b / var_t

  # RIV: relative increase in variance
  riv <- (1 + 1 / M) * var_b / var_w

  # --- Barnard-Rubin adjusted degrees of freedom ---
  if (is.na(v_com)) {
    # v_com is NA: cannot compute df_adj, fmi, or re
    df_adj <- NA_real_
  } else if (is.infinite(v_com) && var_b == 0) {
    # Large-sample with no between-imputation variance
    df_adj <- Inf
  } else {
    # v_old: degrees of freedom from old (Rubin 1987) formula
    if (lambda != 0) {
      v_old <- (M - 1) / lambda^2
    }

    # v_obs: observed-data degrees of freedom (Barnard-Rubin adjustment)
    if (!is.infinite(v_com)) {
      v_obs <- ((v_com + 1) / (v_com + 3)) * v_com * (1 - lambda)
    }

    # Combine v_old and v_obs
    if (lambda != 0) {
      df_adj <- if (is.infinite(v_com)) v_old else (v_old * v_obs) / (v_old + v_obs)
    } else {
      df_adj <- v_obs
    }
  }

  # --- Adjusted FMI (mice convention) ---
  if (is.na(df_adj)) {
    fmi <- NA_real_
  } else {
    fmi <- (riv + 2 / (df_adj + 3)) / (1 + riv)
  }

  # --- Relative efficiency ---
  if (is.na(fmi)) {
    re <- NA_real_
  } else {
    re <- 1 / (1 + fmi / M)
  }

  list(
    var_w  = var_w,
    var_b  = var_b,
    var_t  = var_t,
    lambda = lambda,
    riv    = riv,
    df_adj = df_adj,
    dfcom  = v_com,
    fmi    = fmi,
    re     = re
  )
}
