#' Describe an rbmi Draws Object
#'
#' Extracts structured metadata from an rbmi draws object, including method,
#' formula, sample count, failures, covariance structure, and (for Bayesian
#' methods) MCMC convergence diagnostics. Returns an S3 object with an
#' informative [print()] method.
#'
#' @param draws_obj A `draws` object returned by [rbmi::draws()].
#'
#' @return An S3 object of class `c("describe_draws", "list")` containing:
#' \describe{
#'   \item{method}{Human-readable method name (e.g., "Bayesian (MCMC via Stan)")}
#'   \item{method_class}{Raw class name: "bayes", "approxbayes", or "condmean"}
#'   \item{n_samples}{Total number of samples}
#'   \item{n_failures}{Number of failed samples}
#'   \item{formula}{Deparsed model formula string}
#'   \item{covariance}{Covariance structure (e.g., "us")}
#'   \item{same_cov}{Logical; whether same covariance is used across groups}
#'   \item{condmean_type}{(condmean only) "jackknife" or "bootstrap"}
#'   \item{n_primary}{(condmean only) Always 1}
#'   \item{n_resampled}{(condmean only) Number of resampled draws}
#'   \item{bayes_control}{(bayes only) List with warmup, thin, chains, seed}
#'   \item{mcmc}{(bayes with stanfit only) List with rhat, ess, max_rhat,
#'     min_ess, n_params, converged}
#' }
#'
#' @details
#' For conditional mean methods, the sample count is displayed as "1 + N"
#' matching the rbmi convention where the first sample is the primary (full-data)
#' fit and the remaining N are jackknife or bootstrap resamples.
#'
#' For Bayesian methods, MCMC convergence diagnostics (ESS, Rhat) are extracted
#' from the `stanfit` object when `rstan` is available. The `converged` flag
#' uses the Rhat < 1.1 threshold matching rbmi's own convention.
#'
#' @seealso
#' * [rbmi::draws()] to create draws objects
#' * [rbmi::method_condmean()], [rbmi::method_bayes()], [rbmi::method_approxbayes()]
#'   for method specification
#'
#' @examples
#' \dontrun{
#' library(rbmi)
#' draws_obj <- draws(data, data_ice, vars, method_condmean(type = "jackknife"))
#' desc <- describe_draws(draws_obj)
#' print(desc)
#'
#' # Programmatic access
#' desc$method
#' desc$n_samples
#' desc$n_failures
#' }
#'
#' @export
describe_draws <- function(draws_obj) {
  # Input validation

  if (!inherits(draws_obj, "draws")) {
    cli::cli_abort(
      "{.arg draws_obj} must be a {.cls draws} object from {.fn rbmi::draws}.",
      class = c("rbmiUtils_error_type", "rbmiUtils_error")
    )
  }

  # Extract method info
  method_class <- class(draws_obj$method)[2]

  # Determine human-readable method name
  method_name <- switch(
    method_class,
    "bayes" = "Bayesian (MCMC via Stan)",
    "approxbayes" = "Approximate Bayesian",
    "condmean" = {
      condmean_type <- draws_obj$method$type %||% "unknown"
      paste0("Conditional Mean (", condmean_type, ")")
    },
    "Unknown"
  )

  # Extract core info
  info <- list(
    method = method_name,
    method_class = method_class,
    n_samples = length(draws_obj$samples),
    n_failures = draws_obj$n_failures,
    formula = paste(deparse(draws_obj$formula), collapse = " "),
    covariance = draws_obj$method$covariance,
    same_cov = draws_obj$method$same_cov
  )

  # Condmean-specific fields
  if (method_class == "condmean") {
    info$condmean_type <- draws_obj$method$type %||% "unknown"
    info$n_primary <- 1L
    info$n_resampled <- length(draws_obj$samples) - 1L
  }

  # Bayesian-specific fields
  if (method_class == "bayes") {
    info$bayes_control <- list(
      warmup = draws_obj$method$control$warmup,
      thin = draws_obj$method$control$thin,
      chains = draws_obj$method$control$chains,
      seed = draws_obj$method$control$seed
    )

    # MCMC diagnostics from stanfit (only when rstan available and fit exists)
    if (inherits(draws_obj$fit, "stanfit")) {
      if (requireNamespace("rstan", quietly = TRUE)) {
        stan_summary <- rstan::summary(draws_obj$fit)$summary
        rhat_vals <- stan_summary[, "Rhat"]
        ess_vals <- stan_summary[, "n_eff"]

        # Handle all-NA Rhat case explicitly
        non_na_rhat <- rhat_vals[!is.na(rhat_vals)]
        converged <- if (length(non_na_rhat) == 0) {
          NA
        } else {
          all(non_na_rhat < 1.1)
        }

        info$mcmc <- list(
          rhat = rhat_vals,
          ess = ess_vals,
          max_rhat = max(rhat_vals, na.rm = TRUE),
          min_ess = min(ess_vals, na.rm = TRUE),
          n_params = nrow(stan_summary),
          converged = converged
        )
      } else {
        cli::cli_inform(
          "Install {.pkg rstan} to see MCMC convergence diagnostics.",
          .frequency = "once",
          .frequency_id = "rbmiUtils_rstan_suggest"
        )
      }
    }
  }

  structure(info, class = c("describe_draws", "list"))
}


#' Print Method for describe_draws Objects
#'
#' Displays a formatted summary of a draws description using cli formatting.
#'
#' @param x A `describe_draws` object from [describe_draws()].
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns `x` (for pipe chaining).
#'
#' @examples
#' \dontrun{
#' desc <- describe_draws(draws_obj)
#' print(desc)
#' }
#'
#' @export
print.describe_draws <- function(x, ...) {
  cli::cli_h1("Draws Summary")
  cli::cli_text("{.field Method}: {x$method}")
  cli::cli_text("{.field Formula}: {x$formula}")

  # Display samples with "1 + N" format for condmean
  if (x$method_class == "condmean") {
    n_res <- x$n_resampled
    cm_type <- x$condmean_type
    cli::cli_text(
      "{.field Samples}: 1 + {n_res} (1 primary + {n_res} {cm_type})"
    )
  } else {
    cli::cli_text("{.field Samples}: {x$n_samples}")
  }

  cli::cli_text("{.field Failures}: {x$n_failures}")
  cli::cli_text("{.field Covariance}: {x$covariance}")

  if (!is.null(x$same_cov)) {
    same_label <- if (x$same_cov) "Yes" else "No"
    cli::cli_text("{.field Same covariance across groups}: {same_label}")
  }

  cli::cli_rule()

  # MCMC convergence section
  if (!is.null(x$mcmc)) {
    cli::cli_h2("MCMC Convergence")

    if (is.na(x$mcmc$converged)) {
      cli::cli_alert_warning("Convergence unknown (all Rhat values are NA)")
    } else if (x$mcmc$converged) {
      n_p <- x$mcmc$n_params
      cli::cli_alert_success("All Rhat < 1.1 ({n_p} parameter{?s})")
    } else {
      n_bad <- sum(x$mcmc$rhat >= 1.1, na.rm = TRUE)
      cli::cli_alert_warning("{n_bad} parameter{?s} with Rhat >= 1.1")
    }

    cli::cli_text("{.field Max Rhat}: {round(x$mcmc$max_rhat, 3)}")
    cli::cli_text("{.field Min ESS}: {round(x$mcmc$min_ess, 1)}")
  }

  # Bayes control section (when no MCMC but bayes method)
  if (!is.null(x$bayes_control) && is.null(x$mcmc)) {
    cli::cli_h2("MCMC Settings")
    cli::cli_text("{.field Warmup}: {x$bayes_control$warmup}")
    cli::cli_text("{.field Thin}: {x$bayes_control$thin}")
    cli::cli_text("{.field Chains}: {x$bayes_control$chains}")
    if (!is.null(x$bayes_control$seed)) {
      cli::cli_text("{.field Seed}: {x$bayes_control$seed}")
    }
  }

  invisible(x)
}
