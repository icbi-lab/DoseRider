#' Compute Derivatives and Identify Zero Points
#'
#' This function calculates the first and second derivatives of the predicted
#' values from a mixed model and identifies points where these derivatives
#' are approximately zero.
#'
#' @param smooth_pathway Data frame containing smoothed trend predictions, typically from a dose-response model.
#' @param dose_var The name of the dose variable in `smooth_pathway`.
#' @return A list containing the first derivative, second derivative, and
#'         the zero points for both first and second derivatives.
compute_derivatives <- function(smooth_pathway, dose_var) {
  # Agregate by Dose and fit
  mean_trend <- aggregate(predictions ~ Dose, data = as.data.frame(smooth_pathway), FUN = mean)

  # Perform prediction
  preds <- mean_trend$predictions

  # Calculate derivatives
  first_deriv <- diff(preds) / diff(mean_trend[[dose_var]])
  second_deriv <- diff(first_deriv) / diff(mean_trend[[dose_var]][-length(mean_trend[[dose_var]])])

  # Identify zero points
  tolerance <- 0.01*max(mean_trend[[dose_var]])  # Adjust tolerance as needed
  zero_points_first_deriv <- mean_trend[[dose_var]][which(abs(first_deriv) < tolerance)]
  zero_points_second_deriv <- mean_trend[[dose_var]][-length(mean_trend[[dose_var]])][which(abs(second_deriv) < tolerance)]

  return(list(first_derivative = first_deriv,
       second_derivative = second_deriv,
       zero_points_first_deriv = zero_points_first_deriv,
       zero_points_second_deriv = zero_points_second_deriv))
}


#' Compute Benchmark Dose (BMD) from Smoothed Trend
#'
#' This function calculates the Benchmark Dose (BMD) based on a smoothed trend from model predictions.
#' The BMD is identified as the dose level where the predicted response exceeds a threshold defined
#' as a specified number of standard deviations (`z`) above the control response. It's designed to
#' work with smoothed trend data, typically derived from `lmer` or `glmer` model predictions.
#'
#' @param smooth_pathway Data frame containing smoothed trend predictions, typically from a dose-response model.
#' @param dose_var The name of the dose variable in `smooth_pathway`.
#' @param z A numeric value specifying the number of standard deviations above the control response
#'          to define the target response for BMD. Default is 1.
#' @return A list or a numeric value representing the calculated BMD. If no BMD is found within a
#'         reasonable range, NA is returned.
compute_bmd_from_main_trend <- function(smooth_pathway, dose_var, z = 1) {
  # Agregate by Dose and fit
  mean_trend <- aggregate(predictions ~ Dose, data = as.data.frame(smooth_pathway), FUN = mean)


  # Aggregate predictions across all genes (if necessary)
  # Assuming predictions are already aggregated or represent a main trend
  min_dose <- min(smooth_pathway[[dose_var]])
  # Calculate control response (mean response at the lowest dose)
  control_response <- smooth_pathway[smooth_pathway[dose_var] == min_dose,"predictions"]

  # Calculate target response
  mean_control <- mean(control_response)
  mean_SD <- sd(control_response)  # Standard deviation of the predictions
  target_response <- mean_control + z * mean_SD

  # Compute the absolute differences between each predicted response and the target response
  differences <- abs(mean_trend$predictions - target_response)
  tolerance <- max(smooth_pathway[[dose_var]])*0.1  # Adjust tolerance as needed

  # Get the dose corresponding to the smallest difference
  closest_dose <- mean_trend[differences < tolerance,dose_var]
  closest_dose <- closest_dose[closest_dose >= min_dose]
  if (length(closest_dose) > 0){
  # Create a data frame to store the BMD result
  bmd_result <- list(BMD_zSD = closest_dose)
  } else{
    bmd_result = NA
  }
  return(bmd_result)
}
