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
#' @seealso [tidy_pool_obj()], [cards::as_card()], [cards::check_ard_structure()]
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
  if (!requireNamespace("cards", quietly = TRUE)) {
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
