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
  mean_trend <- aggregate(predictions ~ get(dose_var, smooth_pathway), data = as.data.frame(smooth_pathway), FUN = mean)
  colnames(mean_trend) <- c(dose_var, "predictions")
  # Perform prediction
  preds <- mean_trend$predictions

  # Calculate derivatives
  first_deriv <- diff(preds) / diff(mean_trend[[dose_var]])
  second_deriv <- diff(first_deriv) / diff(mean_trend[[dose_var]][-length(mean_trend[[dose_var]])])

  # Identify zero points
  tolerance <- 0.01*max(preds)  # Adjust tolerance as needed
  zero_points_first_deriv <- mean_trend[[dose_var]][which(abs(first_deriv) < tolerance)]
  zero_points_second_deriv <- mean_trend[[dose_var]][-length(mean_trend[[dose_var]])][which(abs(second_deriv) < tolerance)]

  return(list(zero_points_first_deriv = zero_points_first_deriv,
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
compute_bmd_from_main_trend <- function(smooth_pathway, dose_var, z = 1, center_values = T) {
  # Agregate by Dose and fit
  # If center_values option is enabled, adjust predictions
  if (center_values) {
    smooth_pathway <- smooth_pathway %>%
      group_by(gene) %>%
      mutate(predictions = predictions - mean(predictions, na.rm = TRUE)) %>%
      ungroup()
    smooth_pathway <- as.data.frame(smooth_pathway)
  }

  mean_trend <- aggregate(as.formula(paste0("predictions ~ ",dose_var,"")), data = as.data.frame(smooth_pathway), FUN = mean)

  # Aggregate predictions across all genes (if necessary)
  # Assuming predictions are already aggregated or represent a main trend
  min_dose <- min(smooth_pathway[[dose_var]])
  # Calculate control response (mean response at the lowest dose)
  control_response <- smooth_pathway[smooth_pathway[dose_var] == min_dose,"predictions"]

  # Calculate target response
  mean_control <- mean(control_response, na.rm = T)
  mean_SD <- sd(control_response)  # Standard deviation of the predictions
  threshold_up <- mean_control + z * mean_SD
  threshold_down <- mean_control - z * mean_SD

  # Compute the absolute differences between each predicted response and the target response
  # Find doses where response crosses the threshold
  bmd_doses <- c()
  for (i in 3:nrow(mean_trend)) {
    current_dose <- mean_trend[i, dose_var]
    previous_dose <- mean_trend[i - 1, dose_var]
    current_response <- mean_trend[i, "predictions"]
    previous_response <- mean_trend[i - 1, "predictions"]

    # Check for crossing from below to above the upper threshold
    if ((previous_response < threshold_up && current_response >= threshold_up) ||
        (previous_response > threshold_up && current_response <= threshold_up)) {
      ratio_up <- (threshold_up - previous_response) / (current_response - previous_response)
      interpolated_dose_up <- previous_dose + ratio_up * (current_dose - previous_dose)
      bmd_doses <- c(bmd_doses, interpolated_dose_up)
    }

    # Check for crossing from above to below the lower threshold
    if ((previous_response < threshold_down && current_response >= threshold_down) ||
        (previous_response > threshold_down && current_response <= threshold_down))  {
      ratio_down <- (threshold_down - previous_response) / (current_response - previous_response)
      interpolated_dose_down <- previous_dose + ratio_down * (current_dose - previous_dose)
      bmd_doses <- c(bmd_doses, interpolated_dose_down)
    }

  }

  # Return BMD values
  if (length(bmd_doses) > 0) {
    return(bmd_doses)
  } else {
    return(NA)
  }
}


#' Compute Benchmark Dose (BMD) for a Specific Gene within a Gene Set
#'
#' This function calculates the Benchmark Dose (BMD) for a specific gene within a given gene set based on smoothed trend data from DoseRider analysis.
#' The BMD is identified as the dose level where the predicted response for the specific gene exceeds a threshold defined as a specified number of standard deviations (`z`) above the control response.
#'
#' @param dose_rider_results A list containing the results of DoseRider analysis.
#' @param gene_set_name The name of the gene set.
#' @param gene_name The name of the specific gene within the gene set.
#' @param z A numeric value specifying the number of standard deviations above the control response to define the target response for BMD. Default is 1.
#' @return A numeric value representing the calculated BMD for the specified gene. If no BMD is found, NA is returned.
compute_bmd_for_gene_in_geneset <- function(dose_rider_results, gene_set_name, gene_name, z = 1,dose_var = "Dose") {
  if (!gene_set_name %in% names(dose_rider_results)) {
    stop("Specified gene set not found in dose_rider_results.")
  }

  gene_set_data <- dose_rider_results[[gene_set_name]]$Smooth_Predictions[[1]]

  # Check if gene is in the gene set
  if (!gene_name %in% gene_set_data$gene) {
    stop("Specified gene not found in the gene set.")
  }

  # Extract smooth pathway data for the gene
  smooth_pathway <- subset(gene_set_data, gene == gene_name)

  # Calculate BMD as before
  compute_bmd_from_main_trend(smooth_pathway = smooth_pathway, dose_var = dose_var , z = z)

  return(bmd_result)
}

extract_bmd_for_pathway <- function(dose_rider_results, gene_set_name){
  # Initialize an empty list to store BMD values
  all_bmd_values <- c()

  # Extract the specific gene set results from dose_rider_results
  gene_set_results <- dose_rider_results[[gene_set_name]]

  # Check if cluster-specific results exist
  if (!is.null(gene_set_results$ClusterSpecificResults)) {
    # Iterate through each cluster to extract BMD values
    for (cluster_name in names(gene_set_results$ClusterSpecificResults)) {
      cluster_data <- gene_set_results$ClusterSpecificResults[[cluster_name]]
      all_bmd_values <- c(all_bmd_values,cluster_data$BMD)
    }
  } else {
    # Handle case where no cluster-specific results are available
    all_bmd_values <- NA
  }

  return(all_bmd_values)
}

#' Plot Benchmark Dose (BMD) Density and Peaks
#'
#' This function creates a plot visualizing the density of BMD values and highlights
#' the peaks where the highest density of BMD values are found.
#'
#' @param bmd_range_output A list containing the output from `get_bmd_range` function,
#' which includes x (BMD values), y (density), and bmd (peaks).
#'
#' @return A ggplot object visualizing the density of BMD values with peaks marked.
#'
#' @examples
#' bmd_range_output <- get_bmd_range(dose_rider_results)
#' plot_bmd_density_and_peaks(bmd_range_output)
#'
#' @export
get_bmd_range <- function(dose_rider_results){
  # Load necessary library

  # Find local maxima (peaks) in the density estimate
  find_peaks <- function(density_obj) {
    # The peaks are where the slope changes from positive to negative
    slope <- diff(density_obj$y)
    q55 <- quantile(density_obj$y, 0.55)
    peaks <- which(diff(sign(slope)) == -2) + 1
    peaks <- peaks[density_obj$y[peaks] > q55]
    return( density_obj$x[peaks])
  }


  # Hypothetical list of BMD values for different pathways
  bmd_values <- c()

  for (gene_set_name in names(dose_rider_results)) {
    bmd_values <- c(bmd_values, extract_bmd_for_pathway(dose_rider_results, gene_set_name))

  }


  # Combine all BMD values into a single vector for density estimation
  all_bmds <- unique(bmd_values)
  all_bmds <- all_bmds[!is.na(all_bmds)]

  # Perform kernel density estimation
  bmd_density <- density(all_bmds)
  peaks <- find_peaks(bmd_density)
  res <- list(x = bmd_density$x, y = bmd_density$y, bmd=peaks)
  # Printing the high density range
  return(res)

}
