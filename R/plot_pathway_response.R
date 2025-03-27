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
                                 dose_points = 25, sample_subset_size = 10, interp_knots = NULL,
                                 bdd_knots = NULL, spline_degree = NULL) {
  # Ensure necessary columns are present
  required_cols <- c(dose_col, sample_col)
  if (omic == "rnaseq") required_cols <- c(required_cols, "size_factors")
  if (!all(required_cols %in% names(long_df))) stop("Some required columns not found in data")

  # Prepare parameters for expand.grid
  expand_params <- list(dose = seq(min(long_df[[dose_col]]), max(long_df[[dose_col]]), length.out = dose_points))

  # Add size factor or sample subset
  if (omic == "rnaseq") {
    size_factors <- unique(long_df[["size_factors"]])
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

  # Extract knot information from stored spline_info
  interp_knots <- interp_knots
  bdd_knots <- bdd_knots
  spline_degree <- spline_degree

  # Create new data for prediction
  new_data <- expand.grid(expand_params)
  colnames(new_data) <- c(dose_col, if (omic == "rnaseq") "size_factors" else sample_col, "gene", covariates_cols)

  # Apply the same B-spline transformation to new doses
  new_spline_matrix <- splines::bs(
    new_data[[dose_col]],
    knots = interp_knots,               # Correctly extract knots
    Boundary.knots = bdd_knots,          # Correctly extract boundary knots
    degree = spline_degree,              # Ensure same spline degree
    intercept = FALSE
  )

  # Convert to a DataFrame with correctly named columns
  new_spline_df <- as.data.frame(new_spline_matrix)
  colnames(new_spline_df) <- paste0("CubicSpline_", seq_len(ncol(new_spline_df)))

  # Merge new splines with new_data
  new_data <- cbind(new_data, new_spline_df)

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



#' Plot Pathway Response with Enhanced Visualization at Specific Points
#'
#' This function creates a plot showing the pathway response with enhanced line thickness
#' at specific zero points of first and second derivatives.
#'
#' @param dose_rider_results A list containing the results from the DoseRider analysis.
#' @param gene_set_name The name of the gene set for which to plot the response.
#' @param dose_col The name of the column representing dose information. Default is "Dose".
#' @param center_values Logical, indicating whether to center the prediction values. Default is TRUE.
#' @param legend_position The position of the legend in the plot. Default is "none".
#' @param text_size Numeric, specifying the size of the text in the plot. Default is 4.
#' @param margin_space Numeric, specifying the margin space around the plot. Default is 0.
#' @param model_metrics Logical, indicating whether to include model metrics in the plot. Default is FALSE.
#' @param v_size Numeric, specifying the size of the points where derivatives are zero. Default is 0.5.
#' @param annotate_gene Logical, indicating whether to annotate the gene in the plot. Default is FALSE.
#' @param annotation_text_size Numeric, specifying the size of the annotation text. Default is 5.
#' @param draw_bmd Logical, indicating whether to draw the benchmark dose (BMD) on the plot. Default is TRUE.
#' @param plot_original_data Logical, indicating whether to draw original data, or predict the whole range of doses. Default is FALSE.
#' @param clusterResults Boolean, if TRUE the genes within a gene set will be clustered to find similar expression patterns. Defaults to TRUE.
#' @param scale_values Logical, indicating whether to scale the prediction values. Default is TRUE.
#' @return A ggplot object representing the pathway response plot.
#' @export
plot_pathway_response <- function(dose_rider_results, gene_set_name, dose_col = "Dose",
                                  center_values = TRUE, scale_values = F, legend_position = "none", text_size = 4,
                                  margin_space = 0, model_metrics = FALSE, v_size = 0.5,
                                  annotate_gene = FALSE, annotation_text_size = 5, draw_bmd = TRUE,
                                  plot_original_data = T, clusterResults = TRUE) {

  # Extract the gene set results
  gene_set_results <- extract_gene_set_results(dose_rider_results, gene_set_name, plot_original_data)

  # Process smooth pathway data
  mean_data <- prepare_mean_data(gene_set_results, dose_col = dose_col, center_values = center_values, scale_values = scale_values, clusterResults = clusterResults)

  # Initialize plot with gene-specific trends
  p <- initialize_plot(mean_data, dose_col, model_metrics, gene_set_name)

  # Add cluster-specific trends and BMD lines
  p <- add_cluster_trends_and_bmd(p, gene_set_results, mean_data, dose_col, draw_bmd, v_size, clusterResults, amount = 0.3)

  # Add annotations if specified
  if (annotate_gene) {
    p <- add_gene_annotations(p, mean_data, dose_col, annotation_text_size)
  }

  # Customize the plot theme
  p <- p + theme_dose_rider(legend_position = legend_position, text_size = text_size, margin_space = margin_space)
  p <- p + ggtitle(gene_set_name)

  return(p)
}


#' Plot Top Significant Pathway Responses
#'
#' Creates a combined plot of top significant pathway responses from dose_rider_results.
#'
#' @param dose_rider_results A list containing the results from doseRider analysis.
#' @param top An integer specifying the number of top gene sets to include in the plot. Default is 15.
#' @param order_column A character string specifying the column to use for ordering gene sets in the plot.
#' @param legend_position The position of the legend in the plot.
#' @param text_size The size of the text in the plot.
#' @param margin_space The margin space around the plot.
#' @param scale_values Logical, indicating whether to scale the prediction values. Default is TRUE.

#' @return A combined ggplot object with top significant pathway response plots.
#' @importFrom cowplot plot_grid
#' @export
plot_top_pathway_responses <- function(dose_rider_results, top=6, ncol = 3, order_column = "NegLogPValue", decreasing = F,  dose_col = "Dose",
                                       center_values = TRUE, scale_values = TRUE, legend_position = "none", text_size = 4,
                                       margin_space = 0, model_metrics = FALSE, v_size = 0.5,
                                       annotate_gene = FALSE, annotation_text_size = 5, draw_bmd = TRUE,
                                       plot_original_data = F, clusterResults = F) {

  # Select top gene sets
  top_gene_sets <- get_top_genesets(dose_rider_results = dose_rider_results,
                                    top = top,
                                    decreasing = decreasing,
                                    order_column = order_column )

  # List to store individual plots
  plot_list <- list()

  for (gene_set_name in top_gene_sets) {
    res_path <- dose_rider_results[[gene_set_name]]

    best_model <- res_path$best_model

    if (!is.na(best_model) && best_model != "null") {
      plot_list[[gene_set_name]] <- plot_pathway_response(dose_rider_results,
                                                          gene_set_name = gene_set_name, dose_col = dose_col,
                                                          legend_position = legend_position, text_size=text_size, margin_space = margin_space,
                                                          model_metrics = model_metrics, v_size = v_size,annotate_gene = annotate_gene, annotation_text_size,
                                                          draw_bmd = draw_bmd, plot_original_data = plot_original_data, center_values = center_values, scale_values = scale_values,
                                                          clusterResults = clusterResults)
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
plot_gene_set_random_effects <- function(dose_rider_results, dose_col = "Dose", top = 10, order_column = "NegLogPValue", decreasing = F) {
  all_random_effects <- data.frame(gene = character(),
                                   gene_set = character(),
                                   RandomIntercept = numeric(),
                                   RandomEffect1 = numeric(),
                                   stringsAsFactors = FALSE)
  #Top pathways in function of P-Value
  # Select top gene sets
  top_gene_sets <- get_top_genesets(dose_rider_results = dose_rider_results,
                                    top = top,
                                    decreasing = decreasing,
                                    order_column = order_column )

  for (gene_set_name in top_gene_sets) {
    random_effects <- dose_rider_results[[gene_set_name]]$random_effect

    if (length(random_effects$RandomIntercept) > 1) {
      # Replace with actual way to access the model
      random_effects$gene_set <- gsub("_", " ", gene_set_name)
      random_effects$gene <- rownames(random_effects)

      # Ensure that the columns match
      random_effects <- random_effects[, colnames(all_random_effects)]

      all_random_effects <- rbind(all_random_effects, random_effects)
    }
  }

  # Ridge plot
  all_random_effects$gene_set <- unlist(lapply(all_random_effects$gene_set,function(x){str_wrap(x,width = 35)}))
  p <- ggplot(all_random_effects, aes(x = RandomEffect1, y = gene_set, fill = gene_set)) +
    geom_density_ridges() +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    scale_fill_manual(values = custom_palette) +
    labs(x = expression("Gene-Specific Difference in Dose Effect" ~ (b[1])), y = "") +
    theme_ridges() +
    theme(legend.position = "none") + theme_dose_rider()

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

  random_effects <- dose_rider_results[gene_set_name][[1]]$random_effect

  if (length(random_effects$RandomIntercept) > 1) {
    random_effects$gene_set <- gene_set_name
    random_effects$gene <- rownames(random_effects)

    p <- ggplot(random_effects, aes(x = RandomEffect1, y = RandomIntercept, label = gene)) +
      geom_point(color = "black", fill = "white", shape = 21, size = 3, stroke = 2) +
      geom_label_repel(aes(label = gene), box.padding = 0.35, point.padding = 0.3,
                       size = 3, force = 1) +
      labs(
        x = expression("Gene-Specific Difference in Dose Effect" ~ (b[1])),
        y = expression("Gene-Specific Difference in Baseline Expression" ~ (b[0])),
        title = str_wrap(gene_set_name, 35)
      ) +
      theme_minimal()

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
plot_dotplot_top_pathways <- function(dose_rider_results, top = 10, order_column = "NegLogPValue", pvalue_column = "best_model_pvalue", decreasing = F) {
  # Convert dose_rider_results to dataframe
  dose_rider_df <- as.data.frame.DoseRider(dose_rider_results)

  # Add -log10(p-value) and sort by it
  dose_rider_df$NegLogPValue <- -log10(dose_rider_df[[pvalue_column]])
  # Select top gene sets
  top_gene_sets <- get_top_genesets(dose_rider_results = dose_rider_results,
                                    top = top,
                                    decreasing = decreasing,
                                    order_column = order_column )
  top_pathways_df <- dose_rider_df[dose_rider_df$Geneset %in% top_gene_sets,]
  top_pathways_df$Geneset <- unlist(lapply(top_pathways_df$Geneset, function(x){str_wrap(gsub("_", " ", x),35)}))
  # Create Dot Plot
  dot_plot <- ggplot(top_pathways_df, aes(x = NegLogPValue, y = reorder(Geneset, NegLogPValue), size = Genes, fill = best_model)) +
    geom_point(shape = 21) +  # Add stroke for better visibility
    labs(x = expression("-log"[10] * "(Adjusted non-linear P-Value)"), y = "",
         title = "Top Pathways by NegLogPValue") +  # Title example (adjust as needed)
    scale_fill_manual(
      values = c("non_linear_mixed" = "orange", "non_linear_fixed" = "blue",
                 "linear" = "green", "null" = "red"),
      name = "Best Model"
    ) +
    scale_size_continuous(
      name = "Gene Set Size",
      range = c(3, 8),  # Adjust size range for better visual impact
      breaks = c(5, 50, 100),  # Adjust to more reasonable values based on data
      limits = c(0, 100)  # Set reasonable limits if needed
    ) +
    theme_dose_rider(legend_position = "bottom", fix_ratio = F) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.box = "vertical",  # Ensure legends are stacked vertically
      legend.spacing.y = unit(0.5, 'lines'),  # Add space between legend rows
      plot.title = element_text(hjust = 0.5)  # Center the title
    ) +
    guides(
      size = guide_legend(nrow = 1, order = 1),  # Size legend in first row
      fill = guide_legend(nrow = 1, order = 2)   # Model legend in second row
    )

  return(dot_plot)
}

#' Plot Benchmark Dose (BMD) Density and Peaks
#'
#' This function creates a plot visualizing the density of BMD values and highlights
#' the peaks where the highest density of BMD values are found.
#'
#' @param bmd_range_output A list containing the output from `get_bmd_range` function,
#' which includes x (BMD values), y (density), and bmd (peaks).
#' @param  log_bmd Bool. If true compute the BMD text as 10 exponential
#' @return A ggplot object visualizing the density of BMD values with peaks marked.
#' @import ggplot2
#' @examples
#' bmd_range_output <- get_bmd_range(dose_rider_results)
#' plot_bmd_density_and_peaks(bmd_range_output)
#'
#' @export

plot_bmd_density_and_peaks <- function(bmd_range_output, log_bmd = T) {
  # Convert the list to a dataframe for plotting
  data_to_plot <- data.frame(x = bmd_range_output$x, y = bmd_range_output$y)

  # Create the plot
  p <- ggplot(data_to_plot, aes(x = x, y = y)) +
    geom_line() +
    geom_vline(xintercept = bmd_range_output$bmd, color = "red", linetype = "dashed") +
    labs(x = "", y = "Density", title = "BMD Density and Peaks") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_dose_rider()

  if(log_bmd){
    p <- p + annotate("text", x = bmd_range_output$bmd, y = max(data_to_plot$y) * 0.7,
             label = paste0("BMD: ", round(10**(bmd_range_output$bmd), 3)),
             color = "red", angle = 90, vjust = -0.5) + xlab("log(BMD)")
  } else {
    p <- p + annotate("text", x = bmd_range_output$bmd, y = max(data_to_plot$y) * 0.7,
                      label = paste0("BMD: ", round(bmd_range_output$bmd, 3)),
                      color = "red", angle = 90, vjust = -0.5) + xlab("BMD")
  }

  return(p)
}

#' Plot Bubble Plot of BMD with Confidence Intervals
#'
#' This function creates a bubble plot visualizing the BMD results with confidence intervals.
#' The size of the bubbles represents the median BMD, while the error bars show the confidence intervals.
#'
#' @param bmd_bounds_df A data frame containing the BMD results with columns for lower bound, upper bound, mean BMD, and median BMD.
#' @return A ggplot object representing the bubble plot of BMD with confidence intervals.
#' @import ggplot2
#' @importFrom stringr str_wrap
#' @examples
#' \dontrun{
#'   # Assuming bmd_bounds_df is available
#'   bmd_plot <- plot_bmd_confidence_intervals(bmd_bounds_df, top = 10)
#'   print(bmd_plot)
#' }
#'
#' @export
plot_bmd_confidence_intervals <- function(bmd_bounds_df, top = 10) {
  # Sort by the specified order column
  bmd_bounds_df <- bmd_bounds_df[!is.na(bmd_bounds_df$Median_BMD),]
  # Wrap the pathway names for better readability
  bmd_bounds_df$Geneset <- unlist(lapply(bmd_bounds_df$Geneset, function(x){str_wrap(gsub("_", " ", x), 35)}))

  # Create Bubble Plot
  bmd_plot <- ggplot(bmd_bounds_df, aes(x = Mean_BMD, y = reorder(Geneset, Mean_BMD), fill = Best_Model)) +
    geom_errorbarh(aes(xmin = Lower_Bound, xmax = Upper_Bound), height = 0.2, color = "black", size = 1) +
    geom_point(shape = 21, size = 3) +
    scale_fill_manual(values = c("non_linear_mixed" = "orange","non_linear_fixed" = "blue", "linear" = "green", "null" = "red"), name = "Best Model") +
    labs(x = "Benchmark Dose (BMD)", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_dose_rider(fix_ratio = F)

  return(bmd_plot)
}



#' Plot TCD1 Medians and Confidence Intervals by Pathway
#'
#' This function creates a row-based plot for TCD1 results, where each pathway is a row,
#' the x-axis represents the log dose, and the cluster with the lowest median TCD1
#' and its confidence intervals are visualized.
#'
#' @param tcd_df A data frame containing TCD cluster results, including medians and confidence intervals for each cluster.
#'               Columns must include `Geneset`, `Cluster1_Mean`, `Cluster1_CI_Lower`, `Cluster1_CI_Upper`, etc.
#' @param fontsize Numeric. Font size for the plot text.
#' @return A ggplot object visualizing the TCD1 medians and confidence intervals.
#' @import ggplot2
#' @examples
#' \dontrun{
#'   # Assuming tcd_df is available
#'   tcd_plot <- plot_tcd1_confidence_intervals(tcd_df, fontsize = 12)
#'   print(tcd_plot)
#' }
#' @export
plot_tcd1_confidence_intervals <- function(tcd_df, fontsize = 12) {
  # Extract cluster columns
  cluster_cols <- grep("Cluster\\d+_Mean", names(tcd_df), value = TRUE)
  ci_lower_cols <- grep("Cluster\\d+_CI_Lower", names(tcd_df), value = TRUE)
  ci_upper_cols <- grep("Cluster\\d+_CI_Upper", names(tcd_df), value = TRUE)

  # Reshape data to long format for ggplot
  plot_data <- do.call(rbind, lapply(1:nrow(tcd_df), function(i) {
    geneset <- tcd_df$Geneset[i]
    clusters <- seq_along(cluster_cols)
    data.frame(
      Geneset = geneset,
      Cluster = clusters,
      Median = sapply(clusters, function(j) tcd_df[i, cluster_cols[j]]),
      CI_Lower = sapply(clusters, function(j) tcd_df[i, ci_lower_cols[j]]),
      CI_Upper = sapply(clusters, function(j) tcd_df[i, ci_upper_cols[j]])
    )
  }))

  # Filter out rows with NA medians
  plot_data <- plot_data[!is.na(plot_data$Median),]

  # Select the cluster with the lowest median TCD1 for each pathway
  plot_data <- do.call(rbind, lapply(split(plot_data, plot_data$Geneset), function(df) {
    df[which.min(df$Median), ]
  }))

  # Plot the TCD1 medians and confidence intervals
  tcd1_plot <- ggplot(plot_data, aes(x = Median, y = reorder(Geneset, Median))) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.3, size = 0.8, color = "black") +
    geom_point(size = 3, shape = 21, fill = "orange", color = "black") +
    labs(x = "Log Dose (TCD1)", y = "") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    theme_dose_rider(fix_ratio = FALSE, text_size = fontsize)

  return(tcd1_plot)
}


#' Plot Trend Change Dose (TCD) Density with Zero Points Highlighted
#'
#' This function creates a plot visualizing the density of TCD values and marks
#' the zero points with vertical dashed lines and annotations.
#'
#' @param tcd_range_output A list containing the output from `get_tcd_range` function,
#' which includes `x` (TCD values), `y` (density), and `tcd` (zero points).
#' @param fontsize Numeric. Font size for the plot text.
#'
#' @return A ggplot object visualizing the density of TCD values with zero points marked.
#' @import ggplot2
#' @examples
#' \dontrun{
#' tcd_range_output <- get_tcd_range(dose_rider_results)
#' plot_tcd_density(tcd_range_output, fontsize = 12)
#' }
#'
#' @export
plot_tcd_density <- function(tcd_range_output, fontsize = 12) {
  # Convert the list to a dataframe for plotting
  data_to_plot <- data.frame(x = tcd_range_output$x, y = tcd_range_output$y)
  zero_points <- tcd_range_output$tcd

  # Create the plot
  p <- ggplot(data_to_plot, aes(x = x, y = y)) +
    geom_line() +
    geom_vline(xintercept = zero_points, color = "blue", linetype = "dashed") +
    annotate(
      "text",
      x = zero_points,
      y = max(data_to_plot$y) * 0.7,
      label = paste0("TCD: ", round(10**(zero_points), 2)),
      color = "blue",
      angle = 90,
      vjust = -0.5
    ) +
    labs(x = "TCD", y = "Density", title = "TCD Density") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    theme_dose_rider(text_size = fontsize, fix_ratio = FALSE)

  return(p)
}


