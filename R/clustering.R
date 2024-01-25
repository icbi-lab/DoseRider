#' Transform Smooth Values to Gene vs. Dose Dataframe
#'
#' This function transforms smooth values into a dataframe with genes as rows and doses as columns.
#' Each row (gene) is centered to 0 to normalize the data.
#'
#' @param smooth_values Dataframe containing smooth values from doseRider results.
#' @param dose_col Character string specifying the name of the dose column in `smooth_values`.
#' @return A transformed dataframe with genes as rows, doses as columns, and each row centered to 0.
#' @importFrom tidyr spread
#' @import dplyr
#' @examples
#' \dontrun{
#'   # Assuming 'smooth_values' is available
#'   transformed_data <- transform_smooth_values(smooth_values, dose_col = "Dose")
#' }
transform_smooth_values <- function(smooth_values, dose_col) {
  mean_data <- aggregate(as.formula(paste0("predictions ~ ",dose_col," + gene")), data = smooth_values, FUN = mean)

  # Transforming the data
  transformed_data <- mean_data %>%
    #select(sym(dose_col), gene, predictions) %>%
    spread(key = dose_col, value = predictions)

  rownames(transformed_data) <- transformed_data$gene
  transformed_data$gene <- NULL

  # Center each row (gene) to 0
  transformed_data <- transformed_data - rowMeans(transformed_data, na.rm = TRUE)

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
    km_res <- kmeans(data, centers = k, nstart = 25)
    sil <- silhouette(km_res$cluster, dist(data))
    sil_widths[k-1] <- mean(sil[, "sil_width"])
  }

  # Find the optimal number of clusters
  optimal_clusters <- which.max(sil_widths) + 1
  #optimal_clusters <- 4
  final_km_res <- kmeans(data, centers = optimal_clusters, nstart = 25)

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
