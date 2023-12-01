#' Compute Derivatives and Identify Zero Points
#'
#' This function calculates the first and second derivatives of the predicted
#' values from a mixed model and identifies points where these derivatives
#' are approximately zero. It is designed to work with models fitted using
#' `lmer` or `glmer` from the `lme4` package.
#'
#' @param model A fitted model object, either of class `lmerMod` or `glmerMod`.
#' @param data A dataframe containing the variables used in the model.
#' @param dose_var The name of the dose variable in `data`.
#' @param omic A character string indicating the type of omics data.
#'             Use "rnaseq" for RNA-Seq data; otherwise, specify another type.
#' @param covariates An optional vector of names of covariates in `data` to be
#'                   included in the prediction. Default is NULL.
#' @return A list containing the first derivative, second derivative, and
#'         the zero points for both first and second derivatives.
#' @examples
#' # Assuming `model` is a fitted lmerMod or glmerMod object
#' # `data` is a dataframe with dose, size_factors, and other covariates
#' compute_derivatives(model, data, "dose", "rnaseq", c("cov1", "cov2"))
#' @importFrom lme4 lmer glmer
compute_derivatives <- function(model, data, dose_var, omic, covariates = c()) {
  # Check if model is a linear mixed model or a generalized linear mixed model
  if (!inherits(model, "lmerMod") && !inherits(model, "glmerMod")) {
    stop("Invalid model type")
  }

  # Define the dose range for prediction and derivative calculation
  dose_range <- seq(min(data[[dose_var]]), max(data[[dose_var]]), length.out = 500)

  # Create a new data frame for prediction
  new_data <- data.frame(dose_var = rep(dose_range, each = nrow(data)))
  colnames(new_data) <- dose_var
  # Incorporate covariates into new_data if provided
  if (length(covariates) > 0 && all(covariates %in% names(data))) {
    for (cov in covariates) {
      new_data[[cov]] <- rep(data[[cov]], times = length(dose_range))
    }
  } else if (length(covariates) > 0) {
    stop("Some covariates not found in data")
  }

  if (omic == "rnaseq") {
    # Ensure size_factors are available in 'data'
    if (!"size_factors" %in% names(data)) {
      stop("size_factors column not found in data")
    }
    # Replicate size_factors to match the new_data
    new_data$size_factor <- rep(data$size_factors, each = length(dose_range))
  }

  # Perform prediction
  preds <- predict(model, newdata = new_data, re.form = NA, type = "response")

  # Calculate derivatives
  first_deriv <- diff(preds) / diff(new_data[[dose_var]])
  second_deriv <- diff(first_deriv) / diff(new_data[[dose_var]][-length(new_data[[dose_var]])])

  # Identify zero points
  tolerance <- 1  # Adjust tolerance as needed
  zero_points_first_deriv <- new_data[[dose_var]][which(abs(first_deriv) < tolerance)]
  zero_points_second_deriv <- new_data[[dose_var]][-length(new_data[[dose_var]])][which(abs(second_deriv) < tolerance)]

  return(list(first_derivative = first_deriv,
       second_derivative = second_deriv,
       zero_points_first_deriv = zero_points_first_deriv,
       zero_points_second_deriv = zero_points_second_deriv))
}


#' Compute Benchmark Dose (BMD) from Main Trend
#'
#' This function calculates the Benchmark Dose (BMD) for a given model and dataset.
#' BMD is determined as the dose at which the predicted response is a certain
#' number of standard deviations (defined by `z`) away from the control response.
#' The function is designed to work with models fitted using `lmer` or `glmer`
#' from the `lme4` package and is applicable to linear and generalized linear mixed models.
#'
#' @param model A fitted model object, either of class `lmerMod` or `glmerMod`.
#' @param data A dataframe containing the variables used in the model.
#' @param dose_var The name of the dose variable in `data`.
#' @param omic A character string indicating the type of omics data.
#'             Use "rnaseq" for RNA-Seq data; otherwise, specify another type.
#' @param z The number of standard deviations away from the control response
#'          to define the target response.
#' @return A dataframe containing the Benchmark Dose (BMD) for each specified level of `z`.
#' @examples
#' # Assuming `model` is a fitted lmerMod or glmerMod object
#' # `data` is a dataframe with dose, size_factors, and other covariates
#' compute_bmd_from_main_trend(model, data, "dose", "rnaseq", z = 1)
#' @importFrom lme4 lmer glmer
compute_bmd_from_main_trend <- function(model, data, dose_var, omic, z = 1, covariates = c()) {
  # Check if model is a linear mixed model or a generalized linear mixed model
  if (!inherits(model, "lmerMod") && !inherits(model, "glmerMod")) {
    stop("Invalid model type")
  }

  # Define the dose range for prediction and derivative calculation
  dose_range <- seq(min(data[[dose_var]]), max(data[[dose_var]]), length.out = 1000)

  # Create a new data frame for prediction
  if (omic == "rnaseq") {
    # Ensure size_factors are available in 'data'
    if (!"size_factors" %in% names(data)) {
      stop("size_factors column not found in data")
    }
    new_data <- expand.grid(dose_var = dose_range, size_factor = unique(unlist(data$size_factors)))
    colnames(new_data) <- c(dose_var, "size_factor")
  }

  # Incorporate covariates into new_data if provided
  if (length(covariates) > 0 && all(covariates %in% names(data))) {
    for (cov in covariates) {
      new_data[[cov]] <- rep(data[[cov]], times = length(dose_range))
    }
  } else if (length(covariates) > 0) {
    stop("Some covariates not found in data")
  }

  # Perform prediction
  new_data$preds <- predict(model, newdata = new_data, re.form = NA, type = "response")

  # Aggregate predictions across all genes (if necessary)
  # Assuming predictions are already aggregated or represent a main trend

  # Calculate control response (mean response at the lowest dose)
  control_response <- new_data[new_data[dose_var] == min(dose_range),"preds"]

  # Calculate target response
  mean_control <- mean(control_response)
  mean_SD <- sd(control_response)  # Standard deviation of the predictions
  target_response <- mean_control + z * mean_SD

  # Compute the absolute differences between each predicted response and the target response
  differences <- abs(new_data$preds - target_response)
  tolerance <- max(dose_range)*0.01  # Adjust tolerance as needed

  # Get the dose corresponding to the smallest difference
  closest_dose <- new_data[differences < tolerance,dose_var]
  if (length(closest_dose) > 0){
  # Create a data frame to store the BMD result
  bmd_result <- list(BMD_zSD = closest_dose)
  } else{
    bmd_result = NA
  }
  return(bmd_result)
}
