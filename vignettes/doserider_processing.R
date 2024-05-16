#!/usr/bin/env Rscript
# Load necessary libraries
library(argparser)
library(doseRider)
library(DESeq2)
library(limma)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)

# Function to load data
load_data <- function(expr_data_path, metadata_path) {
  expr_data <- read.csv(expr_data_path, sep = "\t")
  metadata <- read.csv(metadata_path, sep = "\t")
  list(expr_data = expr_data, metadata = metadata)
}

# Function for differential expression analysis

run_diff_expr <- function(data, metadata, omic, apply_filter, dose_col) {
  # Create the SummarizedExperiment object
  se <- SummarizedExperiment(assays = as.matrix(data), colData = metadata)

  if (omic == "rnaseq") {
    se <- estimate_model_parameters(se)
    if (apply_filter) {
      # Setup DESeq2 dataset
      dds <- DESeqDataSet(se, design = as.formula(paste("~", dose_col, "+ I(", dose_col, "^2)")))
      dds <- DESeq(dds)  # Perform DEA
      res <- results(dds, alpha = 0.05)  # Extract results with adjusted p-value
      se <- se[rownames(res[res$pvalue < 0.05, ]), ]  # Filter genes
    }
  } else if (omic == "microarray") {
    if (apply_filter) {
      # Setup linear models for microarray data
      design <- model.matrix(as.formula(paste("~", dose_col, "+ I(", dose_col, "^2)")), data = colData(se))
      fit <- lmFit(assays(se)[[1]], design)
      fit <- eBayes(fit)
      topTable <- topTable(fit, coef = "Dose", adjust.method = "BH", p.value = 0.05)
      se <- se[rownames(topTable), ]  # Filter genes
    }
  }

  return(se)  # Return the SummarizedExperiment object
}

# Function to run doseRider analysis
perform_dose_rider_analysis <- function(data, gmt_file, omic, min_size, max_size,
                                        dose_col, output_path) {
  gmt_data <- filter_gmt_by_size(read_gmt(gmt_file), min_size, max_size)
  dose_rider_results <- DoseRiderParallel(se = data, gmt = gmt_data,
                                          dose_col = "Dose", omic = omic,
                                          minGSsize = min_size, maxGSsize = max_size,
                                          method = "bonferroni", modelType = "LMM", num_cores = 10)
  save(dose_rider_results, file = paste0(output_path, "/doseRider_results.rda"))
  write.table(as.data.frame.DoseRider(dose_rider_results), file = paste0(output_path, "/doseRider_results.tsv"), sep = "\t")
  return(dose_rider_results)
}

# Function to plot results
library(ggplot2)

# Updated function to plot results and save them to files
plot_results <- function(dose_rider_results, save_path, dose_col, top) {
  # Plot 1: Dose Response Heatmap
  p1 <- dose_response_heatmap(dose_rider_results, dose_col = dose_col, top = top, order_column = "best_model_pvalue", decreasing = FALSE)
  plotFile1 <- paste0(save_path, "/plot1.jpeg")
  jpeg(file=plotFile1,width = 600, height = 500, units = "px")
    plot(p1)
  dev.off()

  # Plot 2: Gene Set Random Effects
  p2 <- plot_gene_set_random_effects(dose_rider_results, dose_col = dose_col, order_column = "best_model_pvalue", top = top)
  plotFile2 <- paste0(save_path, "/plot2.jpeg")
  ggsave(plotFile2, plot = p2, width = 10, height = 8, units = "in", dpi = 300)

  # Plot 3: Top Pathway Responses
  p3 <- plot_top_pathway_responses(dose_rider_results, top = 6, ncol = 3, order_column = "best_model_pvalue")
  plotFile3 <- paste0(save_path, "/plot3.jpeg")
  ggsave(plotFile3, plot = p3, width = 12, height = 8, units = "in", dpi = 300)

  # Plot 5: Dotplot Top Pathways
  p5 <- plot_dotplot_top_pathways(dose_rider_results, top = top, order_column = "Genes", decreasing = TRUE)
  plotFile5 <- paste0(save_path, "/plot5.jpeg")
  ggsave(plotFile5, plot = p5, width = 12, height = 8, units = "in", dpi = 300)

  # Data for BMD Range
  data_bmd <- get_bmd_range(dose_rider_results = dose_rider_results)

  # Plot 7: BMD Density and Peaks
  p7 <- plot_bmd_density_and_peaks(data_bmd)
  plotFile7 <- paste0(save_path, "/plot7.jpeg")
  ggsave(plotFile7, plot = p7, width = 12, height = 8, units = "in", dpi = 300)
}


# Main function to execute the workflow
main <- function() {
  # Set up argument parser
  p <- arg_parser("Run doseRider Analysis")
  p <- add_argument(p, "--expr_data", help="Path to expression data file", type="character")
  p <- add_argument(p, "--metadata", help="Path to metadata file", type="character")
  p <- add_argument(p, "--dose", help="Dose Column name in metadata", type="character")
  p <- add_argument(p, "--omic", help="Type of omic data (e.g., 'rnaseq', 'microarray')", type="character")
  p <- add_argument(p, "--gmt_file", help="Path to GMT file", type="character")
  p <- add_argument(p, "--min_size", help="Minimum gene set size", type="integer")
  p <- add_argument(p, "--max_size", help="Maximum gene set size", type="integer")
  p <- add_argument(p, "--output_path", help="Output path for results and plots", type="character")
  p <- add_argument(p, "--filter", help="Whether to apply differential expression filtering (TRUE or FALSE)", default=TRUE, type="logical")
  p <- add_argument(p, "--top", help="Top results pathway to plot", type="integer")

  # Parse arguments
  argv <- parse_args(p)

  # Load data
  data_list <- load_data(argv$expr_data, argv$metadata)

  # Differential expression analysis
  se <- run_diff_expr(data = data_list$expr_data, metadata = data_list$metadata,
                      omic = argv$omic, apply_filter = argv$filter, dose_col = argv$dose)

  # Run doseRider
  dose_rider_results <- perform_dose_rider_analysis(data = se, gmt_file = argv$gmt_file,
                                         omic = argv$omic, min_size = argv$min_size,
                                         max_size = argv$max_size, output_path = argv$output_path,
                                         dose_col = argv$dose)

  # Plot and save results
  plot_results(dose_rider_results = dose_rider_results, save_path = argv$output_path,
               dose_col = argv$dose, top = argv$top)

  # Print completion message
  print("doseRider analysis completed and results saved.")
}

# Execute main function
main()
