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
  tolerance <- 0.001 * max(preds)  # Adjust tolerance as needed
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

#' Compute BMD Bounds using Bootstrapping
#'
#' This function calculates the lower and upper bounds of the Benchmark Dose (BMD)
#' using bootstrapping on the raw data for significant gene sets.
#'
#' @param dose_rider_results A DoseRider object containing results of the analysis.
#' @param dose_col Name of the column representing dose information.
#' @param sample_col Name of the column representing sample information.
#' @param covariates Optional, vector specifying the covariate column(s) in `se`.
#' @param omic Type of omics data, defaults to "rnaseq".
#' @param n_bootstrap The number of bootstrap samples to generate. Default is 1000.
#' @param ci_level Confidence interval level, default is 0.95.
#' @param model_type The type of model to filter by. Options are "all", "linear", "non_linear" (which includes "non_linear_fixed" and "non_linear_mixed").
#' @param filter_type The type of p-value to filter by. Options are "pvalue" for raw p-value and "fdr" for adjusted p-value.
#' @param threshold The p-value threshold for filtering. Defaults to 0.05.
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import progress
#' for Generalized Additive Mixed Models. Defaults to "LMM".
#' @return A list containing the lower and upper BMD bounds for each significant gene set.
#'
#' @examples
#' \dontrun{
#' data("SummarizedExperiment")
#' gmt <- list(geneSet1 = list(genes = c("gene1", "gene2")))
#' results <- DoseRider(se, gmt, "dose", "sample", "covariate", "rnaseq", modelType = "GAMM")
#' bmd_bounds <- compute_bmd_bounds(results, "dose", "sample", n_bootstrap = 1000, ci_level = 0.95)
#' print(bmd_bounds)
#' }
#'
#' @export
compute_bmd_bounds <- function(dose_rider_results, dose_col = "dose", sample_col = "sample",
                               covariates = c(), omic = "rnaseq", n_bootstrap = 100,
                               ci_level = 0.95, filter_type = "fdr", threshold = 0.1,
                               top = 10, model_type = "all") {


  # Initialize df to store BMD bounds
  bmd_bounds_df <- data.frame(
    Geneset = character(),
    Lower_Bound = numeric(),
    Upper_Bound = numeric(),
    Mean_BMD = numeric(),
    Median_BMD = numeric(),
    Best_Model = character(),
    stringsAsFactors = FALSE
  )

  # Get the list of significant gene sets
  significant_gene_sets <- head(filter_DoseRider(dose_rider_results,
                                            filter_type = filter_type,
                                            threshold = threshold,
                                            model_type = model_type), top)

  # Initialize progress bar
  pb <- progress::progress_bar$new(
    format = "  Bootstrapping [:bar] :percent in :elapsed, ETA: :eta",
    total = length(significant_gene_sets) * n_bootstrap, clear = FALSE, width = 60
  )

  # For each significant gene set, perform bootstrap sampling
  for (geneset_name in names(significant_gene_sets)) {

    # Get the original gene set data
    geneset <- significant_gene_sets[[geneset_name]]
    long_df <- geneset$Raw_Values[[1]]
    best_model <- geneset$best_model
    formula <- create_lmm_formula("counts", dose_col, "gene", covariates, best_model, omic)

    # Store BMD values from each bootstrap sample
    bmd_values <- c()

    for (i in seq_len(n_bootstrap)) {
      # Bootstrap sampling: sample with replacement
      # Vectorized approach to create the bootstrap sample indices
      bootstrap_indices <- sample.int(n = nrow(long_df), size = nrow(long_df), replace = TRUE)

      # Create the bootstrap sample using the precomputed indices
      long_df_bootstrap <- long_df[bootstrap_indices, , drop = FALSE]

      # Recompute the DoseRider analysis for the bootstrap sample
      bootstrap_results <- suppressMessages(fit_model_compute_bmd(long_df_bootstrap, formula, omic = omic, clusterResults = FALSE, dose_col = dose_col))

      # Access the BMD results
      bmd_value <- bootstrap_results$AllGenes

      # Extract BMD from the bootstrap results
      if (!is.na(bootstrap_results)) {
        bmd_values <- c(bmd_values, unlist(bmd_value)[1])
      }

      # Update progress bar
      pb$tick()
    }

    # Remove NA values from BMD calculations
    bmd_values <- na.omit(bmd_values)

    # Calculate the lower and upper bounds for the BMD
    lower_bound <- quantile(bmd_values, probs = (1 - ci_level) / 2)
    upper_bound <- quantile(bmd_values, probs = 1 - (1 - ci_level) / 2)
    mean_bmd <- mean(bmd_values)
    median_bmd <- quantile(bmd_values, probs = 0.5)

    # Store the bounds in the results list
    bmd_bounds_df <- rbind(bmd_bounds_df, data.frame(
      Geneset = geneset_name,
      Lower_Bound = lower_bound,
      Upper_Bound = upper_bound,
      Mean_BMD = mean_bmd,
      Median_BMD = median_bmd,
      Best_Model = best_model,
      row.names = c(geneset_name),
      stringsAsFactors = FALSE
    ))
  }

  return(bmd_bounds_df)
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
fit_model_compute_bmd <- function(long_df, formula, omic = "rnaseq", clusterResults = FALSE, dose_col) {

  # Fit the model to the bootstrapped data
  model <- suppressWarnings(fit_lmm(formula, long_df, omic))

  # Compute BMD depending on whether clustering is applied or not
  if (clusterResults) {
    # Cluster genes and compute BMD for each cluster
    max_cluster <- length(unique(long_df$gene)) - 1
    smooth_values <- smooth_pathway_trend(model, long_df, dose_col, "sample", omic, TRUE, dose_points = 50)
    optimal_clusters <- optimal_clusters_silhouette(smooth_values, dose_col, max_clusters = max_cluster)
    clusters <- optimal_clusters$Cluster
    n_cluster <- optimal_clusters$OptimalClusters

    bmd_cluster_results <- list()

    for (j in seq_len(n_cluster)) {
      cluster_genes <- names(clusters[clusters == j])
      smooth_cluster <- smooth_values[smooth_values$gene %in% cluster_genes, ]
      bmd_cluster <- compute_bmd_from_main_trend(smooth_cluster, dose_col, z = 1)
      bmd_cluster_results[[paste("Cluster", j)]] <- bmd_cluster
    }

    return(bmd_cluster_results)

  } else {
    # Compute BMD for the entire gene set
    long_df$predictions <- predict(model, newdata = long_df)
    smooth_values <- smooth_pathway_trend(model, long_df, dose_col, "sample", omic, TRUE, dose_points = 50)
    bmd_all_genes <- compute_bmd_from_main_trend(smooth_values, dose_col, z = 1)
    return(list(AllGenes = bmd_all_genes))
  }
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

