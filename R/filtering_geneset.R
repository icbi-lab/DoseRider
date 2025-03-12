
#' Filter Pathways Based on PCA to Detect Antagonist Patterns
#'
#' This function filters pathways in transcriptomic data by performing Principal Component Analysis (PCA).
#' It determines whether to keep or skip a pathway based on the variance explained by the first principal component (PC1),
#' and by checking for antagonist patterns in the loadings of PC1.
#'
#' @param long_df Dataframe containing the long format transcriptomic data.
#' @param dose_col Character string specifying the name of the dose column in `long_df`. Default is "Dose".
#' @param pca_threshold Numeric value specifying the variance threshold for PC1 to filter pathways. Default is 0.6. Higher values are more restrictive, requiring PC1 to explain a larger portion of variance.
#' @param expression_col Character string specifying the name of the expression column in `long_df`. Default is "counts".
#' @param loading_threshold Numeric value specifying the loading threshold to detect antagonist patterns. Default is 0.5. Lower values are more restrictive, identifying more subtle antagonistic patterns.
#' @param antagonist_threshold Numeric value specifying the threshold for detecting antagonist patterns. Default is 0.5. Lower values are more restrictive, identifying more subtle antagonistic patterns, while higher values are less restrictive.
#' @return Logical value: TRUE if the pathway should be kept (no antagonist pattern detected), FALSE otherwise.
#' @importFrom stats prcomp
#' @examples
#' \dontrun{
#'   # Assuming 'long_df' is available and transform_smooth_values function is defined
#'   should_keep_pathway <- filter_pathway_by_pca(long_df, dose_col = "Dose", pca_threshold = 0.7, expression_col = "counts", loading_threshold = 0.5)
#'   if (should_keep_pathway) {
#'     # Proceed with further analysis
#'   }
#' }
filter_pathway_by_pca <- function(long_df, dose_col = "Dose", pca_threshold = 0.7, expression_col = "counts", loading_threshold = 0.6, antagonist_threshold = 2) {
  # Ensure the transform_smooth_values function is available and works correctly
  pathway_data <- transform_smooth_values(long_df, expression_col = expression_col, dose_col = dose_col)

  # Perform PCA on the transposed data
  pca_results <- prcomp(t(pathway_data), scale. = TRUE)

  # Calculate the variance explained by each principal component
  variance_explained <- (pca_results$sdev^2) / sum(pca_results$sdev^2)

  # Assess the variance explained by the first principal component (PC1)
  pc1_variance <- variance_explained[1]

  # Check if PC1 explains a significant portion of variance
  if (pc1_variance > pca_threshold) {
    # Check loadings for antagonist patterns
    pc1_loadings <- pca_results$rotation[, 1]
    positive_loadings <- sum(pc1_loadings[pc1_loadings > loading_threshold])
    negative_loadings <- sum(abs(pc1_loadings[pc1_loadings < -loading_threshold]))

    # Calculate antagonist score
    antagonist_score <- abs(positive_loadings - negative_loadings) / (positive_loadings + negative_loadings + 1e-10)

    # Determine if there are significant antagonist patterns in loadings
    if (antagonist_score < antagonist_threshold) {
      cat("Skipping pathway: Antagonist pattern detected in PCA loadings\n")
      return(FALSE)
    }
  }

  cat("Keeping pathway: No significant antagonist pattern detected in PCA\n")
  return(TRUE)
}


#' Variance-Based Filtering: Remove Gene Sets with Low Variance
#'
#' @param long_df Dataframe in long format with transcriptomic data.
#' @param expression_col Character, name of the expression column. Default: "counts".
#' @param variance_threshold Numeric, percentile threshold for variance filtering. Default: 0.1 (bottom 10%).
#' @return Logical value: TRUE if the pathway passes filtering, FALSE otherwise.
#' @import dplyr
#' @import tidyr
filter_pathway_by_variance <- function(long_df, expression_col = "counts", variance_threshold = 0.1) {

  # Compute variance for each gene
  gene_variances <- long_df %>%
    group_by(gene) %>%
    summarise(variance = var(counts, na.rm = TRUE)) %>%
    ungroup()

  # Determine the variance cutoff based on the threshold percentile
  var_cutoff <- quantile(gene_variances$variance, variance_threshold)

  # Count how many genes fall below the cutoff
  low_variance_genes <- sum(gene_variances$variance < var_cutoff)

  # If the majority of genes have low variance, exclude the pathway
  if (low_variance_genes / length(gene_variances$gene) > 0.5) {
    cat("Skipping pathway: Majority of genes have low variance\n")
    return(FALSE)
  }

  cat("Keeping pathway: Sufficient variance in gene expression\n")
  return(TRUE)
}




#' Correlation-Based Filtering: Remove Gene Sets with High Antagonistic Patterns
#'
#' @param long_df Dataframe in long format with transcriptomic data.
#' @param expression_col Character, name of the expression column. Default: "counts".
#' @param correlation_threshold Numeric, proportion of negatively correlated genes to filter. Default: 0.5.
#' @return Logical value: TRUE if the pathway passes filtering, FALSE otherwise.
#' @import dplyr
#' @import tidyr
filter_pathway_by_correlation <- function(long_df, expression_col = "counts", correlation_threshold = 0.5) {

  # Pivot the long_df to wide format (genes as rows, samples as columns)
  data_df <- long_df %>%
    select(sample, gene, !!sym(expression_col)) %>%
    pivot_wider(names_from = sample, values_from = !!sym(expression_col), values_fill = NA)

  # Convert to matrix and set row names as gene names
  gene_matrix <- as.matrix(data_df[, -1])
  rownames(gene_matrix) <- data_df$gene

  # Compute pairwise correlations between genes
  correlation_matrix <- cor(t(gene_matrix), use = "pairwise.complete.obs")

  # Count the proportion of negative correlations
  num_genes <- nrow(correlation_matrix)
  total_pairs <- num_genes * (num_genes - 1)  # Total number of gene pairs
  negative_correlation_proportion <- sum(correlation_matrix[lower.tri(correlation_matrix)] < 0.5) / total_pairs

  # If the proportion of negative correlations is too high, exclude the pathway
  if (negative_correlation_proportion > correlation_threshold) {
    cat("Skipping pathway: High proportion of antagonistic gene correlations (", round(negative_correlation_proportion, 3), ")\n")
    return(FALSE)
  }

  cat("Keeping pathway: No significant antagonistic correlations detected\n")
  return(TRUE)
}

#' Generalized Pathway Filtering Function
#'
#' @param long_df Dataframe in long format with transcriptomic data.
#' @param filters List of filtering methods to apply. Options: "pca", "variance", "correlation".
#' @param filter_params Named list containing parameters for each filtering method.
#' @return Logical value: TRUE if the pathway passes all selected filters, FALSE otherwise.
apply_pathway_filters <- function(long_df, filters = c("pca", "variance", "correlation"), filter_params = list()) {

  is_empty_list <- function(lst) {
    return(length(lst) == 0)
  }
  # Default filtering parameters
  default_params <- list(
    pca = list(pca_threshold = 0.7, loading_threshold = 0.6, antagonist_threshold = 2),
    variance = list(variance_threshold = 0.1),
    correlation = list(correlation_threshold = 0.5)
  )

  # Merge default parameters with user-provided parameters
  for (filter in filters) {
    if (filter %in% names(filter_params)) {
      default_params[[filter]] <- modifyList(default_params[[filter]], filter_params[[filter]])
    }
  }

  # Apply PCA filtering if selected
  if ("pca" %in% filters) {
    if (!filter_pathway_by_pca(long_df,
                               pca_threshold = default_params$pca$pca_threshold,
                               loading_threshold = default_params$pca$loading_threshold,
                               antagonist_threshold = default_params$pca$antagonist_threshold)) {
      return(FALSE)
    }
  }

  # Apply variance filtering if selected
  if ("variance" %in% filters) {
    if (!filter_pathway_by_variance(long_df,
                                    variance_threshold = default_params$variance$variance_threshold)) {
      return(FALSE)
    }
  }

  # Apply correlation filtering if selected
  if ("correlation" %in% filters) {
    if (!filter_pathway_by_correlation(long_df,
                                       correlation_threshold = default_params$correlation$correlation_threshold)) {
      return(FALSE)
    }
  }

  # If all filters passed, return TRUE
  return(TRUE)
}


