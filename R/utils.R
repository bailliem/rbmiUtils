#' Get Imputed Data Sets as a data frame
#'
#' This function takes an imputed dataset and a mapping variable to return a dataset
#' with the original IDs mapped back and renamed appropriately.
#'
#' @param impute_obj The imputation object from which the imputed datasets are extracted.
#'
#' @return A data frame with the original subject IDs mapped and renamed.
#' @export
#' @examples
#' \donttest{
#' library(dplyr)
#' library(rbmi)
#' library(rbmiUtils)
#'
#' set.seed(1974)
#' # Load example dataset
#' data("ADEFF")
#'
#' # Prepare data
#' ADEFF <- ADEFF |>
#'   mutate(
#'     TRT = factor(TRT01P, levels = c("Placebo", "Drug A")),
#'     USUBJID = factor(USUBJID),
#'     AVISIT = factor(AVISIT)
#'   )
#'
#' # Define variables for imputation
#' vars <- set_vars(
#'   subjid = "USUBJID",
#'   visit = "AVISIT",
#'   group = "TRT",
#'   outcome = "CHG",
#'   covariates = c("BASE", "STRATA", "REGION")
#' )
#'
#' # Define Bayesian imputation method
#' method <- method_bayes(
#'   n_samples = 100,
#'   control = control_bayes(warmup = 200, thin = 2)
#' )
#'
#' # Generate draws and perform imputation
#' draws_obj <- draws(data = ADEFF, vars = vars, method = method)
#' impute_obj <- impute(draws_obj,
#'   references = c("Placebo" = "Placebo", "Drug A" = "Placebo"))
#'
#' # Extract imputed data with original subject IDs
#' admi <- get_imputed_data(impute_obj)
#' head(admi)
#'}
get_imputed_data <- function(impute_obj) {
  # Check class of impute_obj
  if (!inherits(impute_obj, "imputation")) {
    stop(
      "impute_obj must be of an imputation object outputted from rbmi::impute"
    )
  }

  # Extract the `subjid` variable from the `vars` list
  uid <- impute_obj$data$vars$subjid

  # Extract imputed datasets
  imputed_dfs <- rbmi::extract_imputed_dfs(impute_obj, idmap = TRUE)

  # Extract the ID map from the attributes of the first imputed dataset
  # example code from issue https://github.com/insightsengineering/rbmi/issues/382
  idmap <- attributes(imputed_dfs[[1]])$idmap

  # Convert the imputed data list into a data frame, adding an IMPID variable
  imputed_dfs2 <- imputed_dfs |>
    purrr::map_dfr(~.x, .id = "IMPID")

  # Map original IDs back to the data
  imputed_dfs2$original_id <- idmap[match(imputed_dfs2[[uid]], names(idmap))]

  # Rename to ensure original ID has the correct name
  result <- imputed_dfs2 |>
    dplyr::rename(internal_id = {{ uid }}) |>
    dplyr::rename({{ uid }} := original_id)

  return(result)
}


#' Utility function for Generalized G-computation for Binary Outcomes
#'
#' Wrapper function for targeting a marginal treatment effect
#' using g-computation using the beeca package. Intended for binary endpoints.
#'
#' @param data A data.frame containing the analysis dataset.
#' @param outcome Name of the binary outcome variable (as string).
#' @param treatment Name of the treatment variable (as string).
#' @param covariates Character vector of covariate names to adjust for.
#' @param reference Reference level for the treatment variable (default: "Placebo").
#' @param contrast Type of contrast to compute (default: "diff").
#' @param method Marginal estimation method for variance (default: "Ge").
#' @param type Variance estimator type (default: "HC0").
#' @param ... Additional arguments passed to `beeca::get_marginal_effect()`.
#'
#' @return A named list with treatment effect estimate, standard error, and degrees of freedom (if applicable).
#'
#' @export
#'
#' @examples
#' # Load required packages
#' library(rbmiUtils)
#' library(beeca)      # for get_marginal_effect()
#' library(dplyr)
#' # Load example data
#' data("ADMI")
#' # Ensure correct factor levels
#' ADMI <- ADMI %>%
#'   mutate(
#'     TRT = factor(TRT, levels = c("Placebo", "Drug A")),
#'     STRATA = factor(STRATA),
#'     REGION = factor(REGION)
#'   )
#' # Apply g-computation for binary responder
#' result <- gcomp_binary(
#'   data = ADMI,
#'   outcome = "CRIT1FLN",
#'   treatment = "TRT",
#'   covariates = c("BASE", "STRATA", "REGION"),
#'   reference = "Placebo",
#'   contrast = "diff",
#'   method = "Ge",    # from beeca: GEE robust sandwich estimator
#'   type = "HC0"      # from beeca: heteroskedasticity-consistent SE
#' )
#'
#' # Print results
#' print(result)
#'
gcomp_binary <- function(
  data,
  outcome = "CRIT1FLN",
  treatment = "TRT",
  covariates = c("BASE", "STRATA", "REGION"),
  reference = "Placebo",
  contrast = "diff",
  method = "Ge",
  type = "HC0",
  ...
) {
  # Construct formula
  form <- stats::as.formula(
    paste0(outcome, " ~ ", paste(c(treatment, covariates), collapse = " + "))
  )

  # Fit logistic regression
  model <- stats::glm(form, data = data, family = binomial)

  # Compute marginal treatment effect
  marginal_fit <- beeca::get_marginal_effect(
    model,
    trt = treatment,
    method = method,
    type = type,
    contrast = contrast,
    reference = reference,
    ...
  )

  res <- marginal_fit$marginal_results

  out <- list(
    trt = list(
      est = res[res$STAT == paste0(contrast), "STATVAL"][[1]],
      se = res[res$STAT == paste0(contrast, "_se"), "STATVAL"][[1]],
      df = NA
    )
  )

  return(out)
}


#' Reduce an imputed dataset to imputed values only
#'
#' This function takes a fully imputed dataset and removes rows where the
#' original data was observed. It is useful for storing only the values that
#' were generated during imputation.
#'
#' @param data A data frame containing all imputed datasets. Must include a
#'   variable `IMPID` identifying the imputation number.
#' @param original_data The original dataset before imputation was performed.
#' @param vars Character vector of variables that were imputed.
#' @param keys Optional character vector of columns used to match `data` to
#'   `original_data`. By default the common columns excluding `vars` and `IMPID`.
#'
#' @return A data frame containing only rows with imputed values.
#' @export
reduce_imputed_data <- function(data, original_data, vars, keys = NULL) {
  if (is.null(keys)) {
    keys <- intersect(names(data), names(original_data))
    keys <- setdiff(keys, c(vars, "IMPID"))
  }

  orig_subset <- original_data[, c(keys, vars), drop = FALSE]
  names(orig_subset) <- c(keys, paste0(vars, ".orig"))

  joined <- dplyr::left_join(data, orig_subset, by = keys)

  keep <- apply(
    as.data.frame(lapply(vars, function(v) {
      is.na(joined[[paste0(v, ".orig")]])
    })),
    1,
    any
  )

  out <- joined[keep, , drop = FALSE]
  out <- dplyr::select(out, -dplyr::any_of(paste0(vars, ".orig")))

  return(out)
}


#' Expand a reduced imputed dataset back to full data
#'
#' Given a dataset created by `reduce_imputed_data()`, this function recreates
#' the full set of imputed datasets by merging the imputed values with the
#' original data.
#'
#' @param reduced A reduced imputed dataset produced by
#'   `reduce_imputed_data()`.
#' @param original_data The original dataset before imputation.
#' @param vars Character vector of variables that were imputed.
#' @param keys Optional character vector used to match rows. Defaults to the
#'   common columns excluding `vars` and `IMPID`.
#'
#' @return A data frame containing the full imputed datasets.
#' @export
expand_imputed_data <- function(reduced, original_data, vars, keys = NULL) {
  if (is.null(keys)) {
    keys <- intersect(names(reduced), names(original_data))
    keys <- setdiff(keys, c(vars, "IMPID"))
  }

  imps <- sort(unique(reduced$IMPID))
  out_list <- lapply(imps, function(i) {
    imp_dat <- reduced[reduced$IMPID == i, , drop = FALSE]
    base <- original_data
    base$IMPID <- i
    base <- dplyr::left_join(base, imp_dat, by = c(keys, "IMPID"),
                             suffix = c("", ".imp"))
    for (v in vars) {
      imp_col <- paste0(v, ".imp")
      base[[v]] <- ifelse(is.na(base[[v]]), base[[imp_col]], base[[v]])
      base[[imp_col]] <- NULL
    }
    base
  })

  dplyr::bind_rows(out_list)
}

