#' Generate Smooth Predictions for Pathway Trends
#'
#' This function creates a data frame of smooth predictions based on a provided model. It is useful for
#' visualizing trends in pathway data, especially in the context of dose-response studies. It handles
#' different types of omics data and can incorporate covariates.
#'
#' @param model A model object, typically of class 'lmerMod' or 'glmerMod'.
#' @param long_df A data frame containing the input data with columns corresponding to dose, sample, and other variables.
#' @param dose_col The name of the dose variable in `long_df`.
#' @param sample_col The name of the sample variable in `long_df`.
#' @param omic A character string indicating the type of omics data (e.g., "rnaseq").
#' @param random_effects Logical, indicating if random effects (like gene) should be included.
#' @param covariates_cols Optional; a vector of names of additional covariates in `long_df`.
#'
#' @return A data frame containing the original data along with the predictions from the model.
#'
#' @examples
#' \dontrun{
#' data("mtcars")
#' model <- lmer(mpg ~ wt + (1|cyl), data = mtcars)
#' predictions <- smooth_pathway_trend(model = model,
#'                                    long_df = mtcars,
#'                                    dose_col = "wt",
#'                                    sample_col = "cyl",
#'                                    omic = "rnaseq",
#'                                    random_effects = TRUE,
#'                                    covariates_cols = NULL)
#' }
#'
#' @importFrom lme4 lmer
#' @importFrom stats predict
smooth_pathway_trend <- function(model, long_df, dose_col = "dose", sample_col = "sample", omic = "rnaseq",
                                 random_effects = TRUE, covariates_cols = NULL, gene_subset = NULL,
                                 dose_points = 25, sample_subset_size = 10) {
  # Ensure necessary columns are present
  required_cols <- c(dose_col, sample_col)
  if (omic == "rnaseq") required_cols <- c(required_cols, "size_factor")
  if (!all(required_cols %in% names(long_df))) stop("Some required columns not found in data")

  # Prepare parameters for expand.grid
  expand_params <- list(dose = seq(min(long_df[[dose_col]]), max(long_df[[dose_col]]), length.out = dose_points))

  # Add size factor or sample subset
  if (omic == "rnaseq") {
    size_factors <- unique(long_df[["size_factor"]])
    if (length(size_factors) > sample_subset_size) {
      size_factors <- sample(size_factors, sample_subset_size)
    }
    expand_params$size_factor <- size_factors
  } else {
    samples <- unique(long_df[[sample_col]])
    if (length(samples) > sample_subset_size) {
      samples <- sample(samples, sample_subset_size)
    }
    expand_params$sample <- samples
  }

  # Add gene if random effects are considered
  if (random_effects) {
    expand_params$gene <- unique(long_df$gene)
  } else {
    expand_params$gene <- unique(long_df$gene)[1]

  }
  # Add additional covariates if provided
  if (!is.null(covariates_cols) && all(covariates_cols %in% names(long_df))) {
    for (col in covariates_cols) expand_params[[col]] <- unique(long_df[[col]])
  } else if (!is.null(covariates_cols)) stop("Some covariates not found in data")

  # Create new data for prediction
  new_data <- expand.grid(expand_params)
  colnames(new_data) <- c(dose_col, if (omic == "rnaseq") "size_factor" else sample_col, "gene", covariates_cols)

  # Generate predictions
  if (random_effects){
    predictions <- predict(model, newdata = new_data)
  } else {
    predictions <- predict(model, newdata = new_data, re.form = NA)
  }
  predictions <- cbind(new_data, predictions)

  # Return the data with predictions
  return(predictions)
}


#Function to clos the TCD trend
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
#'
#' @import ggplot2
#' @importFrom stringr str_wrap
#'
#' @export
plot_pathway_response <- function(dose_rider_results, gene_set_name, dose_col = "Dose", center_values = T) {

  # Extract the specific gene set results from dose_rider_results
  gene_set_results <- dose_rider_results[[gene_set_name]]

  # Smooth Predictions for the gene set
  smooth_pathway <- gene_set_results$Smooth_Predictions[[1]]

  # Cluster Assignments for the gene set, converted to a data frame and adding gene names as a column
  ClusterAssignments <- data.frame(Cluster = gene_set_results$ClusterAssignments)
  ClusterAssignments$gene <- rownames(ClusterAssignments)

  # Cluster-specific results and the number of optimal clusters for the gene set
  ClusterSpecificResults <- gene_set_results$ClusterSpecificResults
  n_cluster <- gene_set_results$OptimalClusters

  # If there are smoothed data
  if (!is.null(smooth_pathway) && length(smooth_pathway) > 0) {
    mean_data <- aggregate(as.formula(paste0("predictions ~ ",dose_col," + gene")), data = as.data.frame(smooth_pathway), FUN = mean)
    mean_data <- merge(mean_data, ClusterAssignments, by = "gene")

    # If center_values option is enabled, adjust predictions
    if (center_values) {
      mean_data <- mean_data %>%
        group_by(gene) %>%
        mutate(predictions = predictions - mean(predictions, na.rm = TRUE)) %>%
        ungroup()
      predictions_col <- "predictions_centered"
    } else {
      predictions_col <- "predictions"
    }

    # Now proceed with plotting,

    # Create custom legend labels
    legend_labels <- create_legend_labels(dose_rider_results, gene_set_name)

    # Add gene-specific trends
    p <- ggplot() + geom_line(data = mean_data, aes(x = Dose, y = predictions, group = gene,color = as.factor(Cluster)), linewidth = 0.5)
    #PLot specific data for the clusters
    for (j in c(1:n_cluster)) {
      cluster_data <- ClusterSpecificResults[paste0("Cluster ",j)][[1]]
      cluster_derivate <- c(cluster_data$Derivative$zero_points_first_deriv,cluster_data$Derivative$zero_points_second_deriv)
      cluster_bmd <- cluster_data$BMD
      cluster_trend <- aggregate(as.formula(paste0("predictions ~ ",dose_col)), data = as.data.frame(mean_data[mean_data$Cluster == j,]), FUN = mean)

      #Plot Cluster trend
      #Plot TCD
      cluster_trend <- adjust_trend_visuals(cluster_trend, dose_col, cluster_derivate)
      p <- p + geom_line(data = cluster_trend, aes(x = Dose, y = predictions, linewidth = line_thickness))
      #Plot BMD
      if (sum(!is.na(cluster_bmd))>0){
        for (bmd_dose in cluster_bmd) {

          p <- p + geom_vline(xintercept = bmd_dose,color="#F57C00",linetype="dashed",linewidth=1.5)

        }
      }


    }
    p <- p + scale_linewidth(range = c(0.7, 1.6))
    # Create the base plot with mean trend
    #p <- p +
    #  geom_line(data = mean_trend, aes(x = Dose, y = predictions,linewidth = line_thickness))

    # Finalize the plot settings
    # Finalize the plot settings and add custom legend labels
    p <- p + labs(x = "Dose", y = "Expression", title = str_wrap(gene_set_name, width = 35), caption = paste(legend_labels)) +
      theme_minimal() +
      theme(
        legend.position = "none", # Place the legend at the bottom
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16)
      )

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
#' @param top An integer specifying the number of top gene sets to include in the plot. Default is 15.
#' @param order_column A character string specifying the column to use for ordering gene sets in the plot.
#'
#'
#' @return A combined ggplot object with top significant pathway response plots.
#' @importFrom cowplot plot_grid
#' @export
plot_top_pathway_responses <- function(dose_rider_results, top=15, ncol = 3, order_column = "best_model_pvalue", decreasing = F) {
  # Extract and order gene sets by adjusted cubic p-value
  dose_rider_df <- as.data.frame.DoseRider(dose_rider_results)
  dose_rider_df <- dose_rider_df[order(dose_rider_df[[order_column]], decreasing = decreasing), ]

  # Select top gene sets
  top_gene_sets <- head(dose_rider_df$Geneset, top)

  # List to store individual plots
  plot_list <- list()

  for (gene_set_name in top_gene_sets) {
    res_path <- dose_rider_results[[gene_set_name]]

    best_model <- res_path$best_model

    if (!is.na(best_model) && best_model != "null") {
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
#' @param dose_col A character string specifying the name of the dose column in the raw expression data.
#' @param top An integer specifying the number of top gene sets to include in the plot. Default is 15.
#' @param order_column A character string specifying the column to use for ordering gene sets in the plot
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
plot_gene_set_random_effects <- function(dose_rider_results, dose_col = "Dose", top = 10, order_column = "best_model_pvalue", decreasing = F) {
  all_random_effects <- data.frame(gene = character(),
                                   gene_set = character(),
                                   RandomIntercept = numeric(),
                                   RandomEffect = numeric(),
                                   stringsAsFactors = FALSE)
  #Top pathways in function of P-Value
  dose_rider_df <- as.data.frame(dose_rider_results)
  dose_rider_df <- dose_rider_df[order(dose_rider_df[[order_column]], decreasing = decreasing), ]

  # Extract the top gene set names
  gene_set_names <- dose_rider_df$Geneset[1:top]

  for (gene_set_name in gene_set_names) {
    random_effects <- dose_rider_results[[gene_set_name]]$random_effect

    if (length(random_effects$RandomEffect) > 1){
    # Replace with actual way to access the model
    random_effects$gene_set <- gene_set_name
    random_effects$gene <- rownames(random_effects)
    all_random_effects <- rbind(all_random_effects, random_effects)
    }
  }

  # Ridge plot
  all_random_effects$gene_set <- unlist(lapply(all_random_effects$gene_set,function(x){str_wrap(x,width = 35)}))
  p <- ggplot(all_random_effects, aes(x = RandomEffect1, y = gene_set, fill = gene_set)) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    geom_density_ridges() +
    scale_fill_manual(values = custom_palette) +
    labs(x = "Random Effect", y = "Gene Set") +
    theme_ridges() +
    theme(legend.position = "none")

  return(p)
}


#' Plot Scatter Plot for Relationship Between Random Intercepts and Effects in a Gene Set
#'
#' This function generates a scatter plot to visualize the relationship between
#' random intercepts and random slopes for genes within a specified gene set.
#' It is particularly useful for analyzing the variability and relationship
#' of gene expressions within a gene set modeled.
#'
#' @param dose_rider_results A list containing the results from the DoseRider analysis.
#' @param gene_set_name The name of the gene set for which the plot will be generated.
#'
#' @return A ggplot object representing the scatter plot of random intercepts vs random slopes for the specified gene set.
#' The plot displays the relationship between random intercepts and random slopes for each gene in the gene set.
#'
#' @import ggplot2
#' @import ggrepel
#' @importFrom stringr str_wrap
#'
#' @examples
#' \dontrun{
#'   # Assuming 'dose_rider_results' contains LMM results for gene sets
#'   scatter_plot <- plot_gene_random_effect_relationship(dose_rider_results, "Gene Set Name")
#'   print(scatter_plot)
#' }
#'
#' @export
plot_gene_random_effect_relationship <- function(dose_rider_results, gene_set_name) {
  if (!gene_set_name %in% names(dose_rider_results)) {
    stop("The specified gene set name is not found in the results.")
  }

  random_effects <- dose_rider_results[[gene_set_name]]$random_effect

  if (length(random_effects$RandomEffect) > 1) {
    random_effects$gene_set <- gene_set_name
    random_effects$gene <- rownames(random_effects)

    p <- ggplot(random_effects, aes(x = RandomEffect1, y = RandomIntercept, label = gene)) +
      geom_point(color = "black", fill = "white", shape = 21, size = 3, stroke = 2) +
      geom_label_repel(aes(label = gene), box.padding = 0.35, point.padding = 0.3,
                       size = 3, force = 1) +
      labs(x = "Random Slope", y = "Random Intercept", title = str_wrap(gene_set_name, 35)) +
      theme_bw()

    return(p)
  } else {
    stop("Insufficient data for random effects in the specified gene set.")
  }
}

#' Plot Dot Plot of Top Pathways from DoseRider Results
#'
#' This function creates a dot plot visualizing the top pathways from DoseRider analysis results based on their adjusted cubic p-value.
#' The size of the dots represents the number of genes in each pathway, while the color indicates the best model selected for each pathway.
#'
#' @param dose_rider_results A list containing the results from the DoseRider analysis.
#' @param top The number of top pathways to display in the plot. Default is 10.
#' @param order_column A character string specifying the column to use for ordering gene sets in the plot.
#' @return A ggplot object representing the dot plot of top pathways.
#' @import ggplot2
#' @importFrom stringr str_wrap
#' @examples
#' \dontrun{
#'   # Assuming dose_rider_results is available
#'   dot_plot <- plot_dotplot_top_pathways(dose_rider_results, top = 10)
#'   print(dot_plot)
#' }
#'
#' @export
plot_dotplot_top_pathways <- function(dose_rider_results, top = 10, order_column = "best_model_pvalue", pvalue_column = "best_model_pvalue", decreasing = F) {
  # Convert dose_rider_results to dataframe
  dose_rider_df <- as.data.frame.DoseRider(dose_rider_results)

  # Add -log10(p-value) and sort by it
  dose_rider_df$NegLogPValue <- -log10(dose_rider_df[[pvalue_column]])
  dose_rider_df <- dose_rider_df[order(dose_rider_df[[order_column]], decreasing = decreasing), ]
  dose_rider_df <- dose_rider_df[!is.na(dose_rider_df$NegLogPValue),]
  top_pathways_df <- head(dose_rider_df, top)
  top_pathways_df$Geneset <- unlist(lapply(top_pathways_df$Geneset, function(x){str_wrap(x,35)}))
  # Create Dot Plot
  dot_plot <- ggplot(top_pathways_df, aes(x = NegLogPValue, y = reorder(Geneset, NegLogPValue), size = Genes, fill = best_model)) +
    geom_point(shape = 21) +
    scale_size_continuous(name = "Gene Set Size") +
    labs(x = "-log10(Adjusted non-linear P-Value)", y = "") +
    scale_fill_manual(values = c("non_linear_mixed" = "orange","non_linear_fixed" = "blue", "linear" = "green", "null" = "red"), name = "Best Model") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(dot_plot)
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
#' @import ggplot2
#' @examples
#' bmd_range_output <- get_bmd_range(dose_rider_results)
#' plot_bmd_density_and_peaks(bmd_range_output)
#'
#' @export

plot_bmd_density_and_peaks <- function(bmd_range_output) {
  # Convert the list to a dataframe for plotting
  data_to_plot <- data.frame(x = bmd_range_output$x, y = bmd_range_output$y)

  # Create the plot
  p <- ggplot(data_to_plot, aes(x = x, y = y)) +
    geom_line() +
    geom_vline(xintercept = bmd_range_output$bmd, color = "red", linetype = "dashed") +
    labs(x = "BMD", y = "Density", title = "BMD Density and Peaks") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}


