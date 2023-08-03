#' Plot Smooth Curves for Gene Expressions
#'
#' This function plots smooth curves for gene expressions along with confidence intervals.
#' It optionally centers the values around the mean expression for each gene.
#'
#' @param result A list containing the results from the 'doserider' function.
#'               This list must contain a nested list for each gene set, and each nested list
#'               must have a 'predictions' data frame with columns for 'fit', 'se', and 'gene'.
#' @param gene_set_name The name of the gene set to plot.
#' @param center_values Logical indicating whether to center the expression values for each gene
#'                      around the mean expression across all genes. Defaults to FALSE.
#'
#' @return A ggplot2 object displaying the smooth curves for gene expressions, optionally centered.
#' @import ggplot2
#' @importFrom stats poly
#' @examples
#' result <- doserider(...)  # result from doserider function
#' plot_smooth(result, "Gene Set 1", center_values = TRUE)
#'
#' @export
plot_smooth <- function(result, gene_set_name, dose_col="dose", center_values = FALSE) {
  # Get the predictions
  predictions <- as.data.frame(result[[gene_set_name]]$Smooth_Predictions)
  mean_data <- aggregate(fit ~ gene + dose, data = predictions, FUN = mean)

  # If center_values option is enabled, adjust predictions
  if (center_values) {
    # Calculate mean expression and standard deviation across all genes
    mean_expression <- colMeans(mean_data["fit"])
    sd_expression <- apply(mean_data["fit"], 2, sd)

    # Center and scale expression values for each gene around the mean
    mean_data["fit"] <- sweep(mean_data["fit"], 2, mean_expression, "-")
    mean_data["fit"] <- sweep(mean_data["fit"], 2, sd_expression, "/")
  }
  # Plot the smooth curve with confidence interval region
  p <- ggplot(mean_data, aes_string(x = dose_col, y = "fit", color = "gene")) +
    geom_line(aes_string(group = "gene"), color = "blue") +
    geom_smooth(formula = y ~ poly(x, 3), data = predictions, aes_string(x = dose_col, y = "fit"), se = TRUE, method = "lm", color = "red") +
    labs(x = "Dose",
         y = "Normalized Expression",
         title = paste("Gene Set:", gene_set_name)) +  # Add gene set name to the title
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 12),     # Increase the text size of axis labels
          axis.title = element_text(size = 14),    # Increase the text size of axis titles
          plot.title = element_text(size = 16))

  return(p)
}
