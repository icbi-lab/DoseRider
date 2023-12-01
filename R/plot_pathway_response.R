#' Plot Pathway Response with Enhanced Visualization at Specific Points
#'
#' This function creates a plot showing the pathway response with enhanced line thickness
#' at specific zero points of first and second derivatives.
#'
#' @param dose_rider_results A list containing the results from the DoseRider analysis.
#' @param gene_set_name The name of the gene set for which to plot the response.
#' @param dose_col The name of the column representing dose information.
#'
#' @return A ggplot object representing the pathway response plot.
#' @export
plot_pathway_response <- function(dose_rider_results, gene_set_name, dose_col = "Dose", center_values = FALSE) {
  smooth_pathway <- dose_rider_results[[gene_set_name]]$Smooth_Predictions[[1]]
  smooth_trend <- dose_rider_results[[gene_set_name]]$Smooth_Predictions_Pathway[[1]]
  zero_points <- dose_rider_results[[gene_set_name]]$TCD$zero_points_first_deriv

  if (!is.null(smooth_pathway) && length(smooth_pathway) > 0) {
    mean_data <- aggregate(predictions ~ Dose + gene, data = as.data.frame(smooth_pathway), FUN = mean)
    mean_trend <- aggregate(predictions ~ Dose, data = as.data.frame(smooth_trend), FUN = mean)

    # Assuming 'mean_trend' and 'mean_data' contain your data

    # If center_values option is enabled, adjust predictions
    if (center_values) {
      # Calculate mean and standard deviation for centering and scaling
      mean_expression <- mean(mean_data$predictions)
      sd_expression <- sd(mean_data$predictions)

      # Center and scale expression values for each gene
      mean_data$predictions <- (mean_data$predictions - mean_expression) / sd_expression

      # Apply the same process to the mean trend data
      #mean_expression_trend <- mean(mean_trend$predictions)
      #sd_expression_trend <- sd(mean_trend$predictions)

      mean_trend$predictions <- (mean_trend$predictions - mean_expression) / sd_expression
    }


    #Ajust thicknes
    mean_trend$line_thickness <- 1

    # Define the range around zero points
    delta <- max(mean_data$Dose)*0.011
    for (z in unique(zero_points)) {
      range_around_zero <- c(z - delta, z + delta)

      # Find the closest values to the range boundaries in mean_trend$Dose
      closest_lower <- mean_trend$Dose[which.min(abs(mean_trend$Dose - range_around_zero[1]))]
      closest_upper <- mean_trend$Dose[which.min(abs(mean_trend$Dose - range_around_zero[2]))]

      # Check if closest values are within the original range to avoid out-of-range issues
      if (closest_lower >= min(mean_trend$Dose) && closest_upper <= max(mean_trend$Dose)) {
        mean_trend$line_thickness <- ifelse(mean_trend$Dose >= closest_lower & mean_trend$Dose <= closest_upper, 1.5, 1)
      }
    }

    # Now proceed with plotting, using the line_thickness column for adjusting line thickness


    # Add gene-specific trends
    p <- ggplot() + geom_line(data = mean_data, aes(x = Dose, y = predictions, group = gene), color = "blue", linewidth = 0.5)

    # Create the base plot with mean trend
    p <- p +
      geom_line(data = mean_trend, aes(x = Dose, y = predictions,linewidth = line_thickness), color = "red")

    # Finalize the plot settings
    p <- p + labs(x = "Dose", y = "Expression", title = gene_set_name) +
      theme_minimal() +
      theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(size = 16))

    return(p)
  } else {
    stop("No valid smooth pathway predictions found for the specified gene set.")
  }
}


#' Plot Top Significant Pathway Responses
#'
#' Creates a combined plot of top significant pathway responses from dose_rider_results.
#'
#' @param dose_rider_results A list containing the results from doseRider analysis.
#' @param top The number of top pathways to plot based on significance.
#'
#' @return A combined ggplot object with top significant pathway response plots.
#' @importFrom cowplot plot_grid
#' @export
plot_top_pathway_responses <- function(dose_rider_results, top=15, ncol = 3) {
  # Extract and order gene sets by adjusted cubic p-value
  dose_rider_df <- as.data.frame(dose_rider_results)
  dose_rider_df <- dose_rider_df[order(dose_rider_df$Adjusted_Cubic_P_Value),]

  # Select top gene sets
  top_gene_sets <- head(dose_rider_df$Geneset, top)

  # List to store individual plots
  plot_list <- list()

  for (gene_set_name in top_gene_sets) {
    res_path <- dose_rider_results[[gene_set_name]]

    best_model <- res_path$Best_Model_AICc
    adj_p_val <- res_path$Adjusted_Cubic_P_Value

    if (!is.na(best_model) && !is.na(adj_p_val) && best_model == "cubic" && adj_p_val < 0.01) {
      plot_list[[gene_set_name]] <- plot_pathway_response(dose_rider_results, gene_set_name)
    }
  }

  # Combine the plots
  if (length(plot_list) > 0) {
    combined_plot <- plot_grid(plotlist = plot_list, align = 'v', ncol = ncol)
    return(combined_plot)
  } else {
    stop("No significant pathways found for the top criteria.")
  }
}


#' Plot Ridge Plots for Random Effects in Gene Sets
#'
#' This function generates ridge plots to visualize the distribution of random effects for genes
#' within each gene set. It is particularly useful for analyzing the variability of gene
#' expressions within gene sets modeled using Linear Mixed Models (LMMs).
#'
#' @param dose_rider_results A list containing the results from the DoseRider analysis.
#' Each element of the list should be a sublist representing a gene set, including a fitted LMM object.
#'
#' @return A ggplot object representing the ridge plots of random effects for each gene set.
#' The plot displays the distribution of random effects for each gene set along the y-axis.
#'
#' @import ggplot2
#' @import ggridges
#' @importFrom stringr str_wrap
#'
#' @examples
#' \dontrun{
#'   # Assuming 'dose_rider_results' contains LMM results for each gene set
#'   ridge_plot <- plot_gene_set_random_effects(dose_rider_results)
#'   print(ridge_plot)
#' }
#'
#' @export
plot_gene_set_random_effects <- function(dose_rider_results, dose_col = "Dose", top = 10) {
  all_random_effects <- data.frame(gene = character(),
                                   gene_set = character(),
                                   random_effect = numeric(),
                                   stringsAsFactors = FALSE)
  #Top pathways in function of P-Value
  dose_rider_df <- as.data.frame(dose_rider_results)
  dose_rider_df <- dose_rider_df[order(dose_rider_df$Adjusted_Cubic_P_Value),]

  # Gene set names
  gene_set_names <- dose_rider_df$Geneset[1:top]

  for (gene_set_name in gene_set_names) {
    random_effects <- dose_rider_results[[gene_set_name]]$random_effect

    if (length(random_effects) > 1){
    # Replace with actual way to access the model
    random_effects$gene_set <- gene_set_name
    random_effects$gene <- rownames(random_effects)
    random_effects$`(Intercept)` <- NULL
    colnames(random_effects) <- c("random_effect","gene_set","gene" )
    all_random_effects <- rbind(all_random_effects, random_effects)
    }
  }

  # Ridge plot
  all_random_effects$gene_set <- unlist(lapply(all_random_effects$gene_set,function(x){str_wrap(x,width = 35)}))
  p <- ggplot(all_random_effects, aes(x = random_effect, y = gene_set, fill = gene_set)) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    geom_density_ridges() +
    scale_fill_manual(values = custom_palette) +
    labs(x = "Random Effect", y = "Gene Set") +
    theme_ridges() +
    theme(legend.position = "none")

  return(p)
}




