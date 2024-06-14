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
#'
#' @examples
#' \dontrun{
#' smooth_pathway <- data.frame(Dose = c(0, 1, 2, 3), predictions = c(0.1, 0.5, 0.9, 1.2))
#' derivatives <- compute_derivatives(smooth_pathway, "Dose")
#' print(derivatives)
#' }
#'
#' @export
compute_derivatives <- function(smooth_pathway, dose_var) {
  # Aggregate by Dose and fit
  mean_trend <- aggregate(predictions ~ get(dose_var, smooth_pathway), data = as.data.frame(smooth_pathway), FUN = mean)
  colnames(mean_trend) <- c(dose_var, "predictions")
  # Perform prediction
  preds <- mean_trend$predictions

  # Calculate derivatives
  first_deriv <- diff(preds) / diff(mean_trend[[dose_var]])
  second_deriv <- diff(first_deriv) / diff(mean_trend[[dose_var]][-length(mean_trend[[dose_var]])])

  # Identify zero points
  tolerance <- 0.01 * max(preds)  # Adjust tolerance as needed
  zero_points_first_deriv <- mean_trend[[dose_var]][which(abs(first_deriv) < tolerance)]
  zero_points_second_deriv <- mean_trend[[dose_var]][-length(mean_trend[[dose_var]])][which(abs(second_deriv) < tolerance)]

  return(list(zero_points_first_deriv = zero_points_first_deriv,
              zero_points_second_deriv = zero_points_second_deriv))
}


#' Compute Benchmark Dose (BMD) from Smoothed Trend
#'
#' This function calculates the Benchmark Dose (BMD) based on a smoothed trend from model predictions.
#' The BMD is identified as the dose level where the predicted response exceeds a threshold defined
#' as a specified number of standard deviations (`z`) above the control response.
#'
#' @param smooth_pathway Data frame containing smoothed trend predictions, typically from a dose-response model.
#' @param dose_var The name of the dose variable in `smooth_pathway`.
#' @param z A numeric value specifying the number of standard deviations above the control response
#'          to define the target response for BMD. Default is 1.
#' @param center_values Logical, if TRUE the predictions are centered around their mean. Default is TRUE.
#' @return A numeric value representing the calculated BMD. If no BMD is found within a reasonable range, NA is returned.
#'
#' @examples
#' \dontrun{
#' smooth_pathway <- data.frame(Dose = c(0, 1, 2, 3), predictions = c(0.1, 0.5, 0.9, 1.2))
#' bmd <- compute_bmd_from_main_trend(smooth_pathway, "Dose")
#' print(bmd)
#' }
#'
#' @export
compute_bmd_from_main_trend <- function(smooth_pathway, dose_var, z = 1, center_values = TRUE) {
  # Aggregate by Dose and fit
  if (center_values) {
    smooth_pathway <- smooth_pathway %>%
      group_by(gene) %>%
      mutate(predictions = predictions - mean(predictions, na.rm = TRUE)) %>%
      ungroup()
    smooth_pathway <- as.data.frame(smooth_pathway)
  }

  mean_trend <- aggregate(as.formula(paste0("predictions ~ ", dose_var)), data = as.data.frame(smooth_pathway), FUN = mean)

  # Calculate control response (mean response at the lowest dose)
  min_dose <- min(smooth_pathway[[dose_var]])
  control_response <- smooth_pathway[smooth_pathway[dose_var] == min_dose, "predictions"]
  mean_control <- mean(control_response, na.rm = TRUE)
  mean_SD <- sd(control_response, na.rm = TRUE)
  threshold_up <- mean_control + z * mean_SD
  threshold_down <- mean_control - z * mean_SD

  # Find doses where response crosses the threshold
  bmd_doses <- c()
  for (i in 3:nrow(mean_trend)) {
    current_dose <- mean_trend[i, dose_var]
    previous_dose <- mean_trend[i - 1, dose_var]
    current_response <- mean_trend[i, "predictions"]
    previous_response <- mean_trend[i - 1, "predictions"]

    if ((previous_response < threshold_up && current_response >= threshold_up) ||
        (previous_response > threshold_up && current_response <= threshold_up)) {
      ratio_up <- (threshold_up - previous_response) / (current_response - previous_response)
      interpolated_dose_up <- previous_dose + ratio_up * (current_dose - previous_dose)
      bmd_doses <- c(bmd_doses, interpolated_dose_up)
    }

    if ((previous_response < threshold_down && current_response >= threshold_down) ||
        (previous_response > threshold_down && current_response <= threshold_down)) {
      ratio_down <- (threshold_down - previous_response) / (current_response - previous_response)
      interpolated_dose_down <- previous_dose + ratio_down * (current_dose - previous_dose)
      bmd_doses <- c(bmd_doses, interpolated_dose_down)
    }
  }

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
#' @param dose_var The name of the dose variable in the results. Default is "Dose".
#' @param z A numeric value specifying the number of standard deviations above the control response to define the target response for BMD. Default is 1.
#' @return A numeric value representing the calculated BMD for the specified gene. If no BMD is found, NA is returned.
#'
#' @examples
#' \dontrun{
#' dose_rider_results <- DoseRider(se, gmt, "dose", "sample")
#' bmd <- compute_bmd_for_gene_in_geneset(dose_rider_results, "geneSet1", "gene1")
#' print(bmd)
#' }
#'
#' @export
compute_bmd_for_gene_in_geneset <- function(dose_rider_results, gene_set_name, gene_name, dose_var = "Dose", z = 1) {
  if (!gene_set_name %in% names(dose_rider_results)) {
    stop("Specified gene set not found in dose_rider_results.")
  }

  gene_set_data <- dose_rider_results[[gene_set_name]]$Smooth_Predictions[[1]]

  if (!gene_name %in% gene_set_data$gene) {
    stop("Specified gene not found in the gene set.")
  }

  smooth_pathway <- subset(gene_set_data, gene == gene_name)
  bmd_result <- compute_bmd_from_main_trend(smooth_pathway = smooth_pathway, dose_var = dose_var, z = z)

  return(bmd_result)
}


#' Extract Benchmark Dose (BMD) Values for a Pathway
#'
#' This function extracts the Benchmark Dose (BMD) values for all clusters within a specified pathway from the DoseRider results.
#'
#' @param dose_rider_results A list containing the results of DoseRider analysis.
#' @param gene_set_name The name of the gene set/pathway.
#' @return A numeric vector containing the BMD values for the specified pathway. If no cluster-specific results are available, NA is returned.
#'
#' @examples
#' \dontrun{
#' dose_rider_results <- DoseRider(se, gmt, "dose", "sample")
#' bmd_values <- extract_bmd_for_pathway(dose_rider_results, "geneSet1")
#' print(bmd_values)
#' }
#'
#' @export
extract_bmd_for_pathway <- function(dose_rider_results, gene_set_name) {
  all_bmd_values <- c()

  gene_set_results <- dose_rider_results[[gene_set_name]]

  if (!is.null(gene_set_results$ClusterSpecificResults)) {
    for (cluster_name in names(gene_set_results$ClusterSpecificResults)) {
      cluster_data <- gene_set_results$ClusterSpecificResults[[cluster_name]]
      all_bmd_values <- c(all_bmd_values, cluster_data$BMD)
    }
  } else {
    all_bmd_values <- NA
  }

  return(all_bmd_values)
}


#' Compute BMD Range and Find Peaks
#'
#' This function computes the range of Benchmark Dose (BMD) values from the DoseRider results and identifies peaks in the density of BMD values.
#'
#' @param dose_rider_results A list containing the results of DoseRider analysis.
#' @return A list containing the BMD values (`x`), their density (`y`), and the identified peaks (`bmd`).
#'
#' @examples
#' \dontrun{
#' dose_rider_results <- DoseRider(se, gmt, "dose", "sample")
#' bmd_range <- get_bmd_range(dose_rider_results)
#' print(bmd_range)
#' }
#'
#' @export
get_bmd_range <- function(dose_rider_results) {
  find_peaks <- function(density_obj) {
    slope <- diff(density_obj$y)
    q55 <- quantile(density_obj$y, 0.55)
    peaks <- which(diff(sign(slope)) == -2) + 1
    peaks <- peaks[density_obj$y[peaks] > q55]
    return(density_obj$x[peaks])
  }

  bmd_values <- c()
  for (gene_set_name in names(dose_rider_results)) {
    bmd_values <- c(bmd_values, extract_bmd_for_pathway(dose_rider_results, gene_set_name))
  }

  all_bmds <- unique(bmd_values)
  all_bmds <- all_bmds[!is.na(all_bmds)]

  bmd_density <- density(all_bmds)
  peaks <- find_peaks(bmd_density)
  res <- list(x = bmd_density$x, y = bmd_density$y, bmd = peaks)

  return(res)
}

#' Extract TCD (Compute Derivatives and Identify Zero Points) Values for a Pathway
#'
#' This function extracts the TCD (Threshold Concentration Dose) values for all clusters within a specified pathway from the DoseRider results.
#'
#' @param dose_rider_results A list containing the results of DoseRider analysis.
#' @param gene_set_name The name of the gene set/pathway.
#' @return A numeric vector containing the TCD values for the specified pathway. If no cluster-specific results are available, NA is returned.
#'
#' @examples
#' \dontrun{
#' dose_rider_results <- DoseRider(se, gmt, "dose", "sample")
#' tcd_values <- extract_tcd_for_pathway(dose_rider_results, "geneSet1")
#' print(tcd_values)
#' }
#'
#' @export
extract_tcd_for_pathway <- function(dose_rider_results, gene_set_name) {
  all_tcd_values <- c()

  gene_set_results <- dose_rider_results[gene_set_name]
  gene_set_results <- gene_set_results[[1]]

  if (!is.null(gene_set_results$ClusterSpecificResults)) {
    for (cluster_name in names(gene_set_results$ClusterSpecificResults)) {
      cluster_data <- gene_set_results$ClusterSpecificResults[[cluster_name]]
      all_tcd_values <- c(all_tcd_values, cluster_data$Derivative$zero_points_second_deriv,cluster_data$Derivative$zero_points_first_deriv)
    }
  } else {
    all_tcd_values <- NA
  }

  return(all_tcd_values)
}


#' Compute TCD Range and Find Zero Points
#'
#' This function computes the range of TCD (Threshold Concentration Dose) values from the DoseRider results and identifies zero points in the density of TCD values.
#'
#' @param dose_rider_results A list containing the results of DoseRider analysis.
#' @return A list containing the TCD values (`x`), their density (`y`), and the identified zero points (`tcd`).
#'
#' @examples
#' \dontrun{
#' dose_rider_results <- DoseRider(se, gmt, "dose", "sample")
#' tcd_range <- get_tcd_range(dose_rider_results)
#' print(tcd_range)
#' }
#'
#' @export
get_tcd_range <- function(dose_rider_results) {
  find_zero_points <- function(density_obj) {
    slope <- diff(density_obj$y)
    zero_points <- which(diff(sign(slope)) == -2) + 1
    return(density_obj$x[zero_points])
  }

  tcd_values <- c()
  for (gene_set_name in names(dose_rider_results)) {
    tcd_values <- c(tcd_values, extract_tcd_for_pathway(dose_rider_results, gene_set_name))
  }

  all_tcds <- unique(tcd_values)
  all_tcds <- all_tcds[!is.na(all_tcds)]

  tcd_density <- density(all_tcds)
  zero_points <- find_zero_points(tcd_density)
  res <- list(x = tcd_density$x, y = tcd_density$y, tcd = zero_points)

  return(res)
}

