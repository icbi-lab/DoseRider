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
compute_derivatives <- function(smooth_pathway, dose_var, tolerance = 0.001 ) {
  # Aggregate by Dose and fit
  mean_trend <- aggregate(predictions ~ get(dose_var, smooth_pathway), data = as.data.frame(smooth_pathway), FUN = mean)
  colnames(mean_trend) <- c(dose_var, "predictions")
  # Perform prediction
  preds <- mean_trend$predictions

  # Calculate derivatives
  first_deriv <- diff(preds) / diff(mean_trend[[dose_var]])
  second_deriv <- diff(first_deriv) / diff(mean_trend[[dose_var]][-length(mean_trend[[dose_var]])])

  # Identify zero points
  tolerance <- tolerance * max(preds)  # Adjust tolerance as needed
  zero_points_first_deriv <- mean_trend[[dose_var]][which(abs(first_deriv) < tolerance)]
  zero_points_second_deriv <- mean_trend[[dose_var]][-length(mean_trend[[dose_var]])][which(abs(second_deriv) < tolerance)]

  return(list(zero_points_first_deriv = zero_points_first_deriv,
              zero_points_second_deriv = zero_points_second_deriv))
}


#' Compute Benchmark Dose (BMD) from Smoothed Trend
#'
#' This function calculates the Benchmark Dose (BMD) based on a smoothed trend from model predictions.
#' The BMD is identified as the dose level where the predicted response exceeds a threshold defined
#' as a specified number of standard deviations (z) above the control response.
#'
#' @param smooth_pathway Data frame containing smoothed trend predictions, typically from a dose-response model.
#' @param dose_var The name of the dose variable in smooth_pathway.
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
compute_bmd_from_main_trend <- function(smooth_pathway, dose_var, z = 1, center_values = TRUE, scale_values = F) {
  # Center the predictions if specified
  if (center_values) {
    # Center the predictions by subtracting the mean
    smooth_pathway <- smooth_pathway %>%
      group_by(gene) %>%
      mutate(predictions = predictions - mean(predictions, na.rm = TRUE)) %>%
      ungroup()
  }

  if (scale_values) {
    # Scale the predictions by dividing by the standard deviation
    smooth_pathway <- smooth_pathway %>%
      group_by(gene) %>%
      mutate(predictions = predictions / sd(predictions, na.rm = TRUE)) %>%
      ungroup()
  }

  # Aggregate the data by the dose variable using as.formula
  mean_trend <- aggregate(as.formula(paste0("predictions ~ ", dose_var)), data = smooth_pathway, FUN = mean)

  # Calculate control response (mean response at the lowest dose)
  min_dose <- min(smooth_pathway[[dose_var]])
  control_response <- as.vector(unlist(smooth_pathway[smooth_pathway[[dose_var]] == min_dose, "predictions"]))

  # Compute the mean and standard deviation of the control response
  mean_control <- mean(control_response, na.rm = TRUE)
  mean_SD <- sd(control_response, na.rm = TRUE)

  # Calculate the thresholds
  threshold_up <- mean_control + z * mean_SD
  threshold_down <- mean_control - z * mean_SD

  # Print the calculated thresholds for debugging purposes
  # print(paste("Mean response:", mean_control))
  # print(paste("Threshold Up:", threshold_up))
  # print(paste("Threshold Down:", threshold_down))
  # print(paste("Max value:", max(mean_trend$predictions)))
  # print(paste("Min value:", min(mean_trend$predictions)))


  # Vectorized identification of crossing points
  crossing_up <- which(diff(sign(mean_trend$predictions - threshold_up)) != 0)
  crossing_down <- which(diff(sign(mean_trend$predictions - threshold_down)) != 0)

  # Vectorized calculation of interpolated doses
  interpolated_dose_up <- mean_trend[[dose_var]][crossing_up] +
    (threshold_up - mean_trend$predictions[crossing_up]) /
    (mean_trend$predictions[crossing_up + 1] - mean_trend$predictions[crossing_up]) *
    (mean_trend[[dose_var]][crossing_up + 1] - mean_trend[[dose_var]][crossing_up])

  interpolated_dose_down <- mean_trend[[dose_var]][crossing_down] +
    (threshold_down - mean_trend$predictions[crossing_down]) /
    (mean_trend$predictions[crossing_down + 1] - mean_trend$predictions[crossing_down]) *
    (mean_trend[[dose_var]][crossing_down + 1] - mean_trend[[dose_var]][crossing_down])

  # Combine the results
  bmd_doses <- c(interpolated_dose_up, interpolated_dose_down)

  # Return the BMD doses if any were found, otherwise return NA
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



#' Fit Model and Compute BMD on Bootstrapped Data
#'
#' This function fits a model to a bootstrapped sample and computes the BMD,
#' either for the entire gene set or for clusters within the gene set.
#'
#' @param long_df Data frame containing the raw data for the gene set.
#' @param formula The formula used to fit the model (LMM or GAMM).
#' @param omic Type of omics data, defaults to "rnaseq".
#' @param clusterResults Boolean, if TRUE, compute BMD for clusters; otherwise, for the whole gene set. Default is FALSE.
#' @param dose_col Name of the column representing dose information.
#' @return A list containing the BMD for the bootstrapped sample.
#'
#' @examples
#' \dontrun{
#' long_df_bootstrap <- long_df[sample(nrow(long_df), replace = TRUE), ]
#' bootstrap_results <- fit_model_compute_bmd(long_df_bootstrap, formula, "rnaseq", clusterResults = FALSE, dose_col = "dose")
#' print(bootstrap_results)
#' }
fit_model_compute_bmd <- function(long_df, formula, omic = "rnaseq", clusterResults = FALSE, dose_col, z = 1, keep_min_bmd = TRUE,
                                  interp_knots = NULL, bdd_knots = NULL, spline_degree = NULL) {

  # cat("Entering fit_model_compute_bmd\n")
  # cat("Dose column:", dose_col, "\n")
  # cat("Unique doses in long_df:\n"); print(unique(long_df[[dose_col]]))
  # cat("interp_knots:\n"); print(interp_knots)
  # cat("bdd_knots:\n"); print(bdd_knots)
  # cat("spline_degree:\n"); print(spline_degree)

  # Fit the model to the bootstrapped data
  model <- suppressWarnings(fit_lmm(formula, long_df, omic))

  long_df$predictions <- predict(model, newdata = long_df)

  # Add check before calling smooth_pathway_trend
  if (is.null(spline_degree)) {
    cat("Warning: spline_degree is NULL\n")
  }

  # cat("Calling smooth_pathway_trend with:\n")
  # cat(" - dose_col:", dose_col, "\n")
  # cat(" - spline_degree:", spline_degree, "\n")
  # cat(" - interp_knots:\n"); print(interp_knots)
  # cat(" - bdd_knots:\n"); print(bdd_knots)

  if (is.null(interp_knots) || is.null(bdd_knots) || is.null(spline_degree)) {
    stop("Spline parameters missing: interp_knots, bdd_knots, or spline_degree is NULL.")
  }

  if (length(spline_degree) == 0 || spline_degree < 1) {
    stop("Invalid spline_degree: must be integer >= 1.")
  }


  smooth_values <- smooth_pathway_trend(model, long_df, dose_col, "sample", omic, TRUE, dose_points = 50,
                                        interp_knots = interp_knots, bdd_knots = bdd_knots,
                                        spline_degree = ifelse(is.null(spline_degree), 3, spline_degree))

  bmd_all_genes <- compute_bmd_from_main_trend(smooth_values, dose_col, z = 1)

  derivatives <- compute_derivatives(smooth_values, dose_var = dose_col)
  derivatives <- c(derivatives$zero_points_first_deriv, derivatives$zero_points_second_deriv)
  derivatives <- derivatives[is.finite(derivatives)]

  return(list("bmd" = min(bmd_all_genes), "tcd" = min(derivatives)))
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

#' Extract Minimum TCD1 Value for a Pathway
#'
#' This function extracts the minimum TCD1 (Threshold Concentration Dose, second derivative zero points)
#' value for a specified pathway from the DoseRider results.
#'
#' @param dose_rider_results A list containing the results of DoseRider analysis.
#' @param gene_set_name The name of the gene set/pathway.
#' @return A numeric value representing the minimum TCD1 for the specified pathway.
#'         If no cluster-specific results are available, NA is returned.
#'
#' @examples
#' \dontrun{
#' dose_rider_results <- DoseRider(se, gmt, "dose", "sample")
#' tcd1_value <- extract_tcd_for_pathway(dose_rider_results, "geneSet1")
#' print(tcd1_value)
#' }
#'
#' @export
extract_tcd_for_pathway <- function(dose_rider_results, gene_set_name) {
  gene_set_results <- dose_rider_results[gene_set_name]
  gene_set_results <- gene_set_results[[1]]

  if (!is.null(gene_set_results$ClusterSpecificResults)) {
    all_tcd1_values <- c()
    for (cluster_name in names(gene_set_results$ClusterSpecificResults)) {
      cluster_data <- gene_set_results$ClusterSpecificResults[[cluster_name]]
      tcd1_values <- cluster_data$Derivative$zero_points_second_deriv
      all_tcd1_values <- c(all_tcd1_values, tcd1_values)
    }
    # Return the minimum TCD1 value, or NA if none exist
    return(if (length(all_tcd1_values) > 0) min(all_tcd1_values, na.rm = TRUE) else NA)
  } else {
    return(NA)
  }
}


#' Compute Minimum TCD1 Values and Find Zero Points in Density
#'
#' This function computes the density of minimum TCD1 (Threshold Concentration Dose,
#' second derivative zero points) values from the DoseRider results and identifies
#' zero points in the density.
#'
#' @param dose_rider_results A list containing the results of DoseRider analysis.
#' @return A list containing the TCD density (`x`, `y`) and the identified zero points (`tcd`).
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
  # Helper function to identify zero points in density
  find_zero_points <- function(density_obj) {
    slope <- diff(density_obj$y)
    zero_points <- which(diff(sign(slope)) == -2) + 1
    return(density_obj$x[zero_points])
  }

  # Extract the minimum TCD1 values for all pathways
  tcd1_values <- c()
  for (gene_set_name in names(dose_rider_results)) {
    tcd1_value <- extract_tcd_for_pathway(dose_rider_results, gene_set_name)
    if (!is.na(tcd1_value)) {
      tcd1_values <- c(tcd1_values, tcd1_value)
    }
  }

  # Compute density for unique TCD1 values
  all_tcd1_values <- unique(tcd1_values)
  tcd_density <- density(all_tcd1_values)

  # Identify zero points in the density
  zero_points <- find_zero_points(tcd_density)

  # Return the results
  return(list(x = tcd_density$x, y = tcd_density$y, tcd = zero_points))
}

#' Cluster TCD Values
#'
#' This function performs clustering on TCD values, combining bootstrap-derived TCDs and original TCDs.
#'
#' @param tcd_values A numeric vector of TCD values (from bootstrap samples).
#' @param original_tcds A numeric vector of original TCD values from the main analysis.
#' @param k Optional number of clusters. If NULL, the number of clusters is determined by the number of original TCDs.
#' @param method Clustering method. Options are "kmeans" (default) or "gmm" (Gaussian Mixture Models).
#' @return A list containing cluster assignments, cluster means, standard deviations, and sizes.
cluster_tcd_values <- function(tcd_values, original_tcds, k = NULL, method = "kmeans") {
  # Combine bootstrap-derived and original TCD values
  all_tcds <- c(tcd_values, original_tcds)
  unique_tcds <- unique(all_tcds)  # Ensure no duplicate points

  # Determine the number of clusters
  if (is.null(k)) {
    k <- length(original_tcds)  # Use the number of original TCDs as default
  }

  # Perform clustering
  if (method == "kmeans") {
    clustering <- kmeans(unique_tcds, centers = k)
  } else if (method == "gmm") {
    clustering <- Mclust(unique_tcds, G = k)
  } else {
    stop("Unsupported clustering method. Use 'kmeans' or 'gmm'.")
  }

  # Summarize clusters
  cluster_assignments <- clustering$cluster
  cluster_means <- tapply(unique_tcds, cluster_assignments, mean)
  cluster_sds <- tapply(unique_tcds, cluster_assignments, sd)
  cluster_sizes <- table(cluster_assignments)

  return(list(
    assignments = cluster_assignments,
    means = cluster_means,
    sds = cluster_sds,
    sizes = cluster_sizes
  ))
}


#' Process a Single Gene Set for BMD and TCD Statistics
#'
#' This function processes a single gene set to compute BMD bounds and TCD cluster statistics.
#'
#' @param geneset A single gene set from `dose_rider_results`.
#' @param geneset_name The name of the gene set being processed.
#' @param dose_col Name of the column representing dose information.
#' @param sample_col Name of the column representing sample information.
#' @param covariates Vector specifying the covariate column(s).
#' @param omic Type of omics data, defaults to "rnaseq".
#' @param n_bootstrap The number of bootstrap samples to generate.
#' @param ci_level Confidence interval level, default is 0.95.
#' @param clusterResults Logical, whether to compute cluster-specific results.
#' @param z A parameter passed to the model fitting function.
#' @return A data frame with BMD bounds and TCD cluster statistics for the gene set.
process_single_geneset <- function(geneset, geneset_name, dose_col, sample_col, covariates, omic,
                                   n_bootstrap, ci_level, clusterResults, z, keep_min_bmd = TRUE) {
  cat("Processing geneset:", geneset_name, "\n")

  long_df <- geneset$Raw_Values[[1]]
  best_model <- geneset$best_model
  cluster <- geneset$ClusterAssignments
  cluster_counts <- table(cluster)
  max_cluster <- names(which.max(cluster_counts))

  interp_knots <- geneset$Knots
  bdd_knots <- geneset$Boundary_Knots
  spline_degree <- geneset$Degree

  cat("Spline parameters:\n")
  print(list(interp_knots = interp_knots, bdd_knots = bdd_knots, spline_degree = spline_degree))

  formula <- doseRider:::create_lmm_formula("counts", dose_col, "gene", long_df, covariates, best_model, omic)
  cat("Model formula created:\n"); print(formula)

  if (clusterResults) {
    cat("Clustering enabled. Selecting genes from dominant cluster:", max_cluster, "\n")
    genes_in_max_cluster <- names(cluster[cluster == max_cluster])
    cluster_name <- paste0("Cluster ", max_cluster)
    original_tcds <- unlist(geneset$ClusterSpecificResults[[cluster_name]]$Derivative)
    long_df <- long_df[long_df$gene %in% genes_in_max_cluster, ]
  } else {
    cat("Clustering disabled. Using all genes.\n")
    max_cluster <- NA
    original_tcds <- unlist(geneset$ClusterSpecificResults[["AllGenes"]]$Derivative)
  }

  cat("Beginning bootstrap...\n")
  bootstrap_results <- replicate(n_bootstrap, {
    bootstrap_indices <- sample.int(n = nrow(long_df), size = nrow(long_df) * 0.8, replace = TRUE)
    long_df_bootstrap <- long_df[bootstrap_indices, , drop = FALSE]

    cat("Bootstrap sample - nrow:", nrow(long_df_bootstrap), "\n")
    cat("Unique doses:", unique(long_df_bootstrap[[dose_col]]), "\n")

    result <-
      fit_model_compute_bmd(
        long_df = long_df_bootstrap,
        formula = formula,
        omic = omic,
        clusterResults = clusterResults,
        dose_col = dose_col,
        z = z,
        keep_min_bmd = keep_min_bmd,
        interp_knots = interp_knots,
        bdd_knots = bdd_knots,
        spline_degree = spline_degree
      )

    return(result)
  })

  cat("Finished bootstrap. Extracting BMD and TCD values.\n")

  bmd_values <- unlist(bootstrap_results[1, ])
  tcd_values <- unlist(bootstrap_results[2, ])

  bmd_values <- na.omit(bmd_values)
  bmd_values <- bmd_values[is.finite(bmd_values)]

  tcd_values <- na.omit(tcd_values)
  tcd_values <- tcd_values[is.finite(tcd_values)]


  cat("Computing BMD statistics...\n")
  bmd_df <- compute_bmd_statistics(bmd_values, ci_level = ci_level)

  cat("Computing TCD statistics...\n")
  tcd_df <- compute_bmd_statistics(tcd_values, ci_level = ci_level)

  result <- data.frame(
    Geneset = geneset_name,
    Lower_Bound_BMD = bmd_df$lower_bound,
    Upper_Bound_BMD = bmd_df$upper_bound,
    Mean_BMD = bmd_df$mean_bmd,
    Median_BMD = bmd_df$median_bmd,
    Lower_Bound_TCD = tcd_df$lower_bound,
    Upper_Bound_TCD = tcd_df$upper_bound,
    Mean_TCD = tcd_df$mean_bmd,
    Median_TCD = tcd_df$median_bmd,
    Best_Model = best_model,
    stringsAsFactors = FALSE
  )


  cat("Finished processing geneset:", geneset_name, "\n")
  return(result)
}

#' Compute BMD Bounds in Parallel
#'
#' @inheritParams process_single_geneset
#' @param dose_rider_results A DoseRider object containing results of the analysis.
#' @param num_cores The number of cores to use for parallel processing of gene sets.
#' @return A data frame with BMD bounds and TCD cluster statistics for all gene sets.
#' @export
compute_bmd_bounds_parallel <- function(dose_rider_results, dose_col, sample_col="sample", ci_level = 0.95,
                                        covariates = c(), omic = "rnaseq", n_bootstrap = 1000,
                                        num_cores = 5, clusterResults = FALSE, z = 1, keep_min_bmd = TRUE) {
  # Register the parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = length(dose_rider_results), style = 3)

  # Perform parallel processing
  bmd_bounds_list <- foreach(geneset_idx = seq_along(dose_rider_results), .combine = rbind, .packages = c("lme4", "doseRider", "dplyr")) %dopar% {
    geneset_name <- names(dose_rider_results)[geneset_idx]
    geneset <- dose_rider_results[[geneset_name]]
    if (!is.na(geneset_name)){
      result <- process_single_geneset(geneset, geneset_name, dose_col, sample_col, covariates, omic, n_bootstrap, ci_level, clusterResults, z, keep_min_bmd)
    }
    # Update progress bar
    setTxtProgressBar(pb, geneset_idx)
    result
  }

  # Stop the parallel backend and close the progress bar
  stopCluster(cl)
  close(pb)

  return(bmd_bounds_list)
}


#' Compute BMD Bounds Sequentially
#'
#' @inheritParams process_single_geneset
#' @param dose_rider_results A DoseRider object containing results of the analysis.
#' @return A data frame with BMD bounds and TCD cluster statistics for all gene sets.
#' @export
compute_bmd_bounds_sequential <- function(dose_rider_results, dose_col, sample_col = "sample", ci_level = 0.95,
                                          covariates = c(), omic = "rnaseq", n_bootstrap = 1000,
                                          clusterResults = FALSE, z = 1, keep_min_bmd = TRUE) {
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = length(dose_rider_results), style = 3)

  bmd_bounds_raw <- lapply(seq_along(dose_rider_results), function(geneset_idx) {
    cat("Processing geneset index:", geneset_idx, "\n")
    geneset_name <- names(dose_rider_results)[geneset_idx]
    geneset <- dose_rider_results[[geneset_name]]

    if (!is.na(geneset_name)) {
      cat("Geneset name:", geneset_name, "\n")
      result <- tryCatch({
        process_single_geneset(geneset, geneset_name, dose_col, sample_col,
                               covariates, omic, n_bootstrap, ci_level,
                               clusterResults, z, keep_min_bmd)
      }, error = function(e) {
        cat("Error in geneset", geneset_name, ":", conditionMessage(e), "\n")
        return(NULL)
      })
    } else {
      cat("Skipped NA geneset name at index", geneset_idx, "\n")
      result <- NULL
    }

    setTxtProgressBar(pb, geneset_idx)
    return(result)
  })

  bmd_bounds_list <- do.call(rbind, Filter(Negate(is.null), bmd_bounds_raw))

  # Close the progress bar
  close(pb)

  return(bmd_bounds_list)
}
