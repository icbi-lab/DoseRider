
### Helper Functions ###

#' Extract gene set results from dose_rider_results
extract_gene_set_results <- function(dose_rider_results, gene_set_name, plot_original_data) {
  gene_set_results <- dose_rider_results[[gene_set_name]]
  if (plot_original_data) {
    smooth_pathway <- gene_set_results$Raw_Values[[1]]
  } else {
    smooth_pathway <- gene_set_results$Smooth_Predictions[[1]]
  }
  return(gene_set_results)
}

#########################################################################

#' Prepare mean data and perform centering/standardization if required
prepare_mean_data <- function(gene_set_results, dose_col, center_values, clusterResults) {
  # Cluster Assignments
  ClusterAssignments <- data.frame(Cluster = gene_set_results$ClusterAssignments)
  ClusterAssignments$gene <- rownames(ClusterAssignments)

  if (!clusterResults) {
    ClusterAssignments$Cluster <- "AllGenes"
  }

  # Smooth pathway predictions
  smooth_pathway <- gene_set_results$Smooth_Predictions[[1]]

  # Standardize predictions
  if (center_values) {
    smooth_pathway <- smooth_pathway %>%
      group_by(gene) %>%
      mutate(predictions = (predictions - mean(predictions, na.rm = TRUE)) / sd(predictions, na.rm = TRUE)) %>%
      ungroup()
  }

  # Aggregate data
  mean_data <- aggregate(as.formula(paste0("predictions ~ ", dose_col, " + gene")),
                         data = as.data.frame(smooth_pathway),
                         FUN = median)
  mean_data <- merge(mean_data, ClusterAssignments, by = "gene")

  return(mean_data)
}

###########################################################################################

#' Initialize the ggplot object with gene-specific trends
initialize_plot <- function(mean_data, dose_col, model_metrics, gene_set_name) {
  #Custom palette
  custom_palette <- c(custom_palette, c("#F8766D"))
  names(custom_palette) <- c(as.character(1:(length(custom_palette) - 1)), "AllGenes")

  p <- ggplot() +
    geom_line(data = mean_data, aes(x = .data[[dose_col]], y = predictions, group = gene, color = as.factor(Cluster)), linewidth = 0.5) +
    scale_color_manual(values = custom_palette)

  # Add model metrics if enabled
  if (model_metrics) {
    legend_labels <- create_legend_labels(dose_rider_results, gene_set_name)
    p <- p + labs(caption = paste(legend_labels))
  }

  return(p)
}

########################################################################################
#' Add cluster trends and BMD lines to the plot
add_cluster_trends_and_bmd <- function(p, gene_set_results, mean_data, dose_col, draw_bmd, v_size, clusterResults) {
  ClusterSpecificResults <- gene_set_results$ClusterSpecificResults
  n_cluster <- gene_set_results$OptimalClusters

  for (j in seq_len(n_cluster)) {
    if (!clusterResults) {
      cluster_data <- ClusterSpecificResults["AllGenes"][[1]]
      j <- "AllGenes"
    } else {
      cluster_data <- ClusterSpecificResults[paste0("Cluster ", j)][[1]]
    }

    # Derivatives and BMD
    cluster_derivate <- c(cluster_data$Derivative$zero_points_first_deriv, cluster_data$Derivative$zero_points_second_deriv)
    cluster_bmd <- cluster_data$BMD

    # Plot cluster trends
    cluster_trend <- aggregate(as.formula(paste0("predictions ~ ", dose_col)),
                               data = as.data.frame(mean_data[mean_data$Cluster == j, ]),
                               FUN = mean)
    p <- p + geom_line(data = cluster_trend, aes(x = .data[[dose_col]], y = predictions), linewidth = 1)

    # Plot BMD lines
    if (draw_bmd && sum(!is.na(cluster_bmd)) > 0) {
      for (bmd_dose in cluster_bmd) {
        p <- p + geom_vline(xintercept = bmd_dose, color = "#F57C00", linetype = "dashed", linewidth = v_size)
      }
    }
  }

  return(p)
}

#####################################################################################
#' Add gene annotations to the plot
add_gene_annotations <- function(p, mean_data, dose_col, annotation_text_size) {
  sampled_mean_data <- mean_data %>%
    group_by(gene) %>%
    slice_sample(n = 1) %>%
    ungroup()

  p <- p + geom_text_repel(data = sampled_mean_data, aes(x = .data[[dose_col]], y = predictions, label = gene),
                           size = annotation_text_size, nudge_y = 0.2, direction = "y", segment.color = 'grey50',
                           show.legend = FALSE)

  return(p)
}
######################################################################################

#Function to clos the TCD trend
##############################
adjust_trend_visuals <- function(cluster_trend, dose_col, zero_points, t = 0.0001, line_thickness = 0.3, new_line_thickness = 0.4) {
  # Initialize columns if they don't exist
  if (!"line_thickness" %in% names(cluster_trend)) {
    cluster_trend$line_thickness <- line_thickness # Default line thickness
  }
  if (!"color" %in% names(cluster_trend)) {
    cluster_trend$color <- "#D32F2F" # Default color
  }

  # Define the delta for zero point range calculation
  delta <- max(cluster_trend[[dose_col]]) * t

  for (z in unique(zero_points)) {
    range_around_zero <- c(z - delta, z + delta)

    # Find the closest values within the range
    closest_lower <- cluster_trend[[dose_col]][which.min(abs(cluster_trend[[dose_col]] - range_around_zero[1]))]
    closest_upper <- cluster_trend[[dose_col]][which.min(abs(cluster_trend[[dose_col]] - range_around_zero[2]))]

    # Adjust line thickness and color for doses within the zero point range
    within_range <- cluster_trend[[dose_col]] >= closest_lower & cluster_trend[[dose_col]] <= closest_upper
    cluster_trend$line_thickness[within_range] <- new_line_thickness
    cluster_trend$color[within_range] <- "#7B1FA2"
  }

  return(cluster_trend)
}
##################################################################################
### Custom theme

theme_dose_rider <- function(legend_position = "none", text_size = 5, margin_space = 1, fix_ratio = TRUE) {
  # Base theme
  custom_theme <- theme_minimal() +
    theme(
      text = element_text(size = text_size),
      plot.caption = element_text(size = text_size - 2),
      plot.title = element_text(size = text_size, hjust = 0.5),
      axis.title = element_text(size = text_size),
      axis.text = element_text(size = text_size),
      legend.title = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      legend.position = legend_position,
      plot.margin = margin(t = margin_space,  # Top margin
                           r = margin_space,  # Right margin
                           b = margin_space,  # Bottom margin
                           l = margin_space)
    )

  # Conditionally set aspect ratio
  if (fix_ratio) {
    custom_theme <- custom_theme + theme(aspect.ratio = 1)
  }

  return(custom_theme)
}

### Custom palette

custom_palette <- c(
  "#388E3C", "#FBC02D", "#F57C00", "#D32F2F", "#7B1FA2", "#303F9F",
  "#0288D1", "#00796B", "#689F38", "#AFB42B", "#FFA000", "#E64A19",
  "#5D4037", "#616161", "#455A64", "#C2185B", "#512DA8", "#1976D2",
  "#0097A7", "#00796B", "#8D6E63", "#78909C", "#6D4C41", "#546E7A",
  "#BF360C", "#3E2723", "#0D47A1", "#1B5E20", "#33691E", "#827717",
  "#F57F17", "#FF6F00", "#E65100", "#BF360C", "#3E2723", "#263238",
  "#212121", "#DD2C00", "#FF3D00", "#FF6E40", "#FF9E80", "#6A1B9A",
  "#AB47BC", "#BA68C8", "#CE93D8"
)




