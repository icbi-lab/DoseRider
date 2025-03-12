#' Transform Smooth Values to Gene vs. Dose Dataframe
#'
#' This function transforms smooth values into a dataframe with genes as rows and doses as columns.
#' Each row (gene) is centered and scaled (standardized) to normalize the data.
#'
#' @param smooth_values Dataframe containing smooth values from doseRider results.
#' @param dose_col Character string specifying the name of the dose column in `smooth_values`.
#' @param expression_col Character string specifying the name of the expression column in `smooth_values`.
#' @param center_values Logical, whether to center the gene values (default is TRUE).
#' @param scale_values Logical, whether to scale the gene values (default is TRUE).
#' @return A transformed dataframe with genes as rows, doses as columns, and each row centered and scaled.
#' @importFrom tidyr spread
#' @import dplyr
#' @examples
#' \dontrun{
#'   # Assuming 'smooth_values' is available
#'   transformed_data <- transform_smooth_values(smooth_values, dose_col = "Dose")
#' }
transform_smooth_values <- function(smooth_values, dose_col = "Dose", expression_col = "predictions",
                                    center_values = TRUE, scale_values = TRUE) {

  # Step 1: Aggregate the data by dose and gene, applying mean function
  mean_data <- aggregate(as.formula(paste0(expression_col," ~ ",dose_col," + gene")),
                         data = smooth_values,
                         FUN = mean)

  # Step 2: Spread the data to make it wide (genes as rows, doses as columns)
  transformed_data <- mean_data %>%
    spread(key = !!sym(dose_col), value = !!sym(expression_col))

  # Set the rownames to gene names
  rownames(transformed_data) <- transformed_data$gene
  transformed_data$gene <- NULL  # Remove the gene column as it's now rownames

  # Step 3: Centering and scaling the gene rows if specified
  if (center_values) {
    # Center each row (gene) by subtracting the row mean
    transformed_data <- transformed_data - rowMeans(transformed_data, na.rm = TRUE)
  }

  if (scale_values) {
    # Scale each row (gene) by dividing by the row's standard deviation
    transformed_data <- transformed_data / apply(transformed_data, 1, sd, na.rm = TRUE)
  }

  return(transformed_data)
}


#' Determine Optimal Number of Clusters Using Silhouette Method
#'
#' This function determines the optimal number of clusters for dose-response data using the silhouette method.
#' It returns the optimal number of clusters and the cluster assignments for each gene.
#'
#' @param data Dataframe containing dose-response data.
#' @param dose_col Character string specifying the name of the dose column in `data`.
#' @param max_clusters Maximum number of clusters to consider. Default is 10.
#' @return A list containing the optimal number of clusters and a dataframe of cluster assignments.
#' @importFrom cluster silhouette
#' @importFrom factoextra fviz_nbclust
#' @examples
#' \dontrun{
#'   # Assuming 'data' is available
#'   cluster_results <- optimal_clusters_silhouette(data, dose_col = "Dose")
#' }
optimal_clusters_silhouette <- function(data, dose_col, max_clusters = 10) {
  # Data
  data <- transform_smooth_values(data, dose_col = dose_col)

  sil_widths <- numeric(max_clusters - 1)

  for (k in 2:max_clusters) {
    set.seed(123)  # for reproducibility
    km_res <- kmeans(data, centers = k, nstart = max_clusters)
    sil <- silhouette(km_res$cluster, dist(data))
    sil_widths[k-1] <- mean(sil[, "sil_width"])
  }

  # Find the optimal number of clusters
  optimal_clusters <- which.max(sil_widths) + 1
  #optimal_clusters <- 4
  final_km_res <- kmeans(data, centers = optimal_clusters, nstart = max_clusters)

  # Creating a list of genes and their corresponding cluster assignments
  #cluster_assignments <- data.frame(gene = rownames(data), Cluster = final_km_res$cluster)

  # Plotting Silhouette width for each number of clusters
  # sil_plot <- fviz_nbclust(data, FUN = kmeans, method = "silhouette", k.max = max_clusters) +
  #   geom_vline(xintercept = optimal_clusters, linetype = 2) +
  #   labs(title = "Silhouette Method for Optimal Number of Clusters",
  #        subtitle = paste("Optimal number of clusters:", optimal_clusters)) +
  #   theme_minimal()

  return(list(OptimalClusters = optimal_clusters, gene = rownames(data), Cluster = final_km_res$cluster))
}




#' Cluster Original TCD Values
#'
#' This function clusters the original TCD values and determines the number of clusters to use for bootstrap data.
#' If the number of unique TCD values is less than 3, it assigns the number of clusters as the length of unique TCDs
#' and returns a named vector for clustering.
#'
#' @param original_tcds A numeric vector of original TCD values.
#' @param method Clustering method. Options are "kmeans" (default) or "gmm" (Gaussian Mixture Models).
#' @return A list containing the number of clusters (`k`) and the clustering results.
cluster_original_tcds <- function(original_tcds, method = "kmeans") {
  unique_tcds <- unique(original_tcds)  # Ensure no duplicate points

  # # Handle cases where the number of unique TCDs is less than 3
  # if (length(unique_tcds) < 4) {
    k <- length(unique_tcds)  # Number of clusters equals the number of unique TCDs
    clustering <- setNames(seq_along(unique_tcds), unique_tcds)  # Create a named vector with cluster IDs
  #   return(list(
  #     k = k,
  #     clustering = clustering
  #   ))
  # }
  #
  # # Perform clustering if the number of unique TCDs is 3 or more
  # if (method == "kmeans") {
  #   k <- min(length(unique_tcds), 3)  # Use up to 3 clusters
  #   clustering <- kmeans(unique_tcds, centers = k, nstart = 25)
  #
  # } else if (method == "gmm") {
  #   clustering <- Mclust(unique_tcds)
  #   k <- clustering$G  # Number of clusters in GMM
  #
  # } else {
  #   stop("Unsupported clustering method. Use 'kmeans' or 'gmm'.")
  # }

  return(list(
    k = k,
    clustering = clustering
  ))
}



#' Cluster Bootstrap TCD Values and Return Data Frame
#'
#' This function clusters bootstrap-derived TCD values using the number of clusters determined from the original TCDs,
#' and computes cluster-level statistics (means, standard deviations, sizes, and confidence intervals).
#'
#' @param tcd_values A numeric vector of bootstrap-derived TCD values.
#' @param k The number of clusters to use, determined from the original TCDs.
#' @param method Clustering method. Options are "kmeans" (default) or "gmm" (Gaussian Mixture Models).
#' @param ci_level Confidence interval level (default is 0.95).
#' @return A data frame summarizing cluster-level statistics.
cluster_bootstrap_tcds <- function(tcd_values, k, method = "kmeans", ci_level = 0.95) {
  unique_tcds <- unique(tcd_values)  # Ensure no duplicate points

  # Handle case when unique_tcds is empty
  if (length(unique_tcds) == 0) {
    warning("No unique TCD values found. Returning an empty data frame.")
    return(data.frame(
      Cluster_ID = integer(0),
      Cluster_Mean = numeric(0),
      Cluster_SD = numeric(0),
      Cluster_Size = integer(0),
      CI_Lower = numeric(0),
      CI_Upper = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  # Handle case when unique_tcds has fewer than 3 values
  if (length(unique_tcds) < 3) {
    warning("Fewer than 3 unique TCD values. Treating each value as its own cluster.")
    cluster_summary_df <- data.frame(
      Cluster_ID = seq_along(unique_tcds),
      Cluster_Mean = unique_tcds,
      Cluster_SD = rep(NA, length(unique_tcds)),
      Cluster_Size = rep(1, length(unique_tcds)),
      CI_Lower = unique_tcds,
      CI_Upper = unique_tcds,
      stringsAsFactors = FALSE
    )
    return(cluster_summary_df)
  }

  # Adjust k if it exceeds the number of unique values or is less than 1
  if (k > length(unique_tcds) || k < 1) {
    warning(paste("The number of clusters (k =", k,
                  ") is invalid. Adjusting k to", length(unique_tcds), "."))
    k <- max(1, length(unique_tcds))  # Ensure k is at least 1
  }

  if (method == "kmeans") {
    clustering <- kmeans(unique_tcds, centers = k, nstart = 25)
    cluster_assignments <- clustering$cluster
  } else if (method == "gmm") {
    clustering <- Mclust(unique_tcds, G = k)
    cluster_assignments <- clustering$classification
  } else {
    stop("Unsupported clustering method. Use 'kmeans' or 'gmm'.")
  }

  # Summarize clusters
  cluster_means <- tapply(unique_tcds, cluster_assignments, mean)
  cluster_sds <- tapply(unique_tcds, cluster_assignments, sd, na.rm = TRUE)
  cluster_sizes <- table(cluster_assignments)

  # Compute confidence intervals for each cluster
  cluster_ci <- lapply(seq_along(cluster_means), function(cluster_id) {
    cluster_values <- unique_tcds[cluster_assignments == cluster_id]
    lower_bound <- quantile(cluster_values, probs = (1 - ci_level) / 2, na.rm = TRUE)
    upper_bound <- quantile(cluster_values, probs = 1 - (1 - ci_level) / 2, na.rm = TRUE)
    return(c(Lower = lower_bound, Upper = upper_bound))
  })

  # Combine CI results into a matrix
  cluster_ci_matrix <- do.call(rbind, cluster_ci)
  colnames(cluster_ci_matrix) <- c("Lower", "Upper")

  # Create a data frame with cluster-level statistics
  cluster_summary_df <- data.frame(
    Cluster_ID = seq_along(cluster_means),
    Cluster_Mean = as.numeric(cluster_means),
    Cluster_SD = as.numeric(cluster_sds),
    Cluster_Size = as.numeric(cluster_sizes),
    CI_Lower = cluster_ci_matrix[, "Lower"],
    CI_Upper = cluster_ci_matrix[, "Upper"],
    stringsAsFactors = FALSE
  )

  return(cluster_summary_df)
}


#' Compute BMD Statistics
#'
#' This function calculates key statistics for BMD values, including lower and upper bounds,
#' mean, and median, based on a given confidence interval level.
#'
#' @param bmd_values A numeric vector of BMD values.
#' @param ci_level Confidence interval level (default is 0.95).
#' @return A named list containing `lower_bound`, `upper_bound`, `mean_bmd`, and `median_bmd`.
compute_bmd_statistics <- function(bmd_values, ci_level = 0.95) {
  if (length(bmd_values) > 0) {
    lower_bound <- quantile(bmd_values, probs = (1 - ci_level) / 2, na.rm = TRUE)
    upper_bound <- quantile(bmd_values, probs = 1 - (1 - ci_level) / 2, na.rm = TRUE)
    mean_bmd <- mean(bmd_values, na.rm = TRUE)
    median_bmd <- median(bmd_values, na.rm = TRUE)
  } else {
    lower_bound <- NA
    upper_bound <- NA
    mean_bmd <- NA
    median_bmd <- NA
  }

  return(list(
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    mean_bmd = mean_bmd,
    median_bmd = median_bmd
  ))
}

