# Title: Using doseRider for Studying Non-Linear Dose-Response in BPA
# Author: Pablo Monfort-Lanzas
# Date: 2023-08-01

# Load necessary libraries
suppressMessages({
  suppressWarnings({
    library(doseRider)
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    library(edgeR)
    library(doParallel)
  })
})

setwd("/home/monfortl/doseRider")

# Load preprocessed data
load("data/PRJNA869442.rda")

# Define the doses for BPAF compound

# Filter data for BPAF
compound_doses <- list(
  "BPAF" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5)
)
doses <- c(0,compound_doses["BPAF"][[1]])

# Example filtering for BPA
bpaf_data <- PRJNA869442[,(colData(PRJNA869442)$Chemical %in% c("Cells, no treatment","BPAF")) &
                           (colData(PRJNA869442)$Dose %in% doses)]
bpaf_data$sample <- colnames(bpaf_data)

# Prepare the data for DESeq2 analysis
colData(bpaf_data)$Dose <- unlist(colData(bpaf_data)$Dose)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = assay(bpaf_data),
  colData   = colData(bpaf_data),
  design    = ~ Dose * Dose
)

# Run DESeq2 analysis
dds <- DESeq(dds, parallel = F, quiet = TRUE)
res <- results(dds)

# Filter results based on criteria
filter_res <- as.data.frame(res[(res$baseMean > 11) & (res$pvalue < 0.1),])
bpaf_data <- estimate_model_parameters(bpaf_data)

# Filter out low express genes
filter_bpaf_data <- bpaf_data[rownames(bpaf_data) %in% rownames(filter_res),]

# Save filtered metadata and expression data
metadata <- colData(filter_bpaf_data)
expression_data <- assay(filter_bpaf_data)

write.table(metadata, "metadata.tsv", sep="\t", quote=FALSE)
write.table(expression_data, "expression_data.tsv", sep="\t", quote=FALSE)

# Specify the path to your GMT files
gmt_h_path <- "external/h.all.v2023.2.Hs.symbols.gmt"
gmt_perturbation_path <- "/home/monfortl/doseRider/external/c2.cgp.v2023.2.Hs.symbols.gmt"

# Save doseRider results and BMD bounds
geneset_name <- "C2_CGP"  # Replace with actual gene set name
minGSsize = 10
maxGSsize = 1000
geneset_size <- paste0(minGSsize,"_",maxGSsize)  # Calculate gene set size
filter_type <- "fdr"  # Example filter type
threshold <- 0.001  # Example threshold

# Read and filter GMT files
gmt <- filter_gmt_by_size(read_gmt(gmt_perturbation_path), minGSsize, maxGSsize)

# Run doseRider analysis
dose_rider_results <- DoseRiderParallel(
  se = filter_bpaf_data,
  gmt = gmt,
  dose_col = "Dose",
  omic = "rnaseq",
  minGSsize = minGSsize,
  maxGSsize = maxGSsize,
  method = "bonferroni",
  covariates = c(),
  modelType = "LMM",
  num_cores = 10,
  clusterResults = F,
  FilterPathway = F
)

# Convert results to data frame
res_df <- as.data.frame.DoseRider(dose_rider_results)
table(res_df$best_model)

# Filter doseRider results based on FDR
dose_rider_results_filter <- filter_DoseRider(
  dose_rider_results,
  model_type = "all",
  filter_type = filter_type,
  threshold = threshold
)

# Compute BMD bounds
bmd_bounds_df <- doseRider::compute_bmd_bounds_parallel(
  dose_rider_results = dose_rider_results_filter,
  dose_col = "Dose",
  sample_col = "sample",
  covariates = c(),
  omic = "rnaseq",
  n_bootstrap = 100,
  num_cores = 20
)



# Define filenames with metadata
dose_rider_results_filename <- paste0("doseRider_results_", geneset_name, "_size", geneset_size, "_cluster_false.rda")
dose_rider_results_filter_filename <- paste0("doseRider_results_filter_", geneset_name,"_",filter_type, "_thresh", threshold, "_cluster_false.rda")
bmd_bounds_filename <- paste0("bmd_bounds_", geneset_name, "_size", geneset_size, ".rda")

# Save results with metadata in filenames
save(dose_rider_results, file = dose_rider_results_filename)
save(dose_rider_results_filter, file = dose_rider_results_filter_filename)
save(bmd_bounds_df, file = bmd_bounds_filename)


# Define plotting parameters
top <- 20
save_path <- "plots/"
dir.create(save_path, showWarnings = FALSE)
PLOT_WIDTH <- 2500
PLOT_HEIGHT <- 2500

# Generate and save plots
p1 <- dose_response_heatmap(dose_rider_results_filter, dose_col = "Dose", top = top, order_column = "best_model_pvalue", decreasing = FALSE)
jpeg(file=paste0(save_path,"plot1.jpeg"), width = PLOT_WIDTH - 1500, height = PLOT_HEIGHT - 1500, units = "px")
plot(p1, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", padding = unit(c(1, 1, 1, 1), "cm"))
dev.off()

p2 <- plot_gene_set_random_effects(dose_rider_results_filter, dose_col = "Dose", order_column = "best_model_pvalue", top = top) + theme_dose_rider(fix_ratio = F)
ggsave(paste0(save_path,"plot2.jpeg"), plot = p2, width = PLOT_WIDTH, height = PLOT_HEIGHT, units = "px", dpi = 600)

p3 <- plot_top_pathway_responses(dose_rider_results_filter, top = 6, ncol = 2, order_column = "best_model_pvalue", text_size = 5)
ggsave(paste0(save_path,"plot3.jpeg"), plot = p3, width = PLOT_WIDTH, height = PLOT_HEIGHT + 1000, units = "px", dpi = 600)

p4 <- plot_gene_random_effect_relationship(dose_rider_results_filter, "BENPORATH_ES_CORE_NINE_CORRELATED")
ggsave(paste0(save_path,"plot4.jpeg"), plot = p4, width = PLOT_WIDTH + 1800, height = PLOT_HEIGHT + 600, units = "px", dpi = 600)

p5 <- plot_dotplot_top_pathways(dose_rider_results_filter, top = top, order_column = "Genes", decreasing = TRUE)
ggsave(paste0(save_path,"plot5.jpeg"), plot = p5, width = PLOT_WIDTH + 500, height = PLOT_HEIGHT + 500, units = "px", dpi = 600)

p6 <- create_gene_heatmap(dose_rider_results_filter, dose_col = "Dose", gene_set_name = "MENSE_HYPOXIA_UP")
jpeg(file=paste0(save_path,"plot6.jpeg"), width = PLOT_WIDTH - 1500, height = PLOT_HEIGHT - 1500, units = "px")
plot(p6, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", padding = unit(c(1, 1, 1, 1), "cm"))
dev.off()

data_bmd <- get_bmd_range(dose_rider_results = dose_rider_results_filter)
p7 <- plot_bmd_density_and_peaks(data_bmd)
ggsave(paste0(save_path,"plot7.jpeg"), plot = p7, width = PLOT_WIDTH, height = PLOT_HEIGHT, units = "px", dpi = 600)

p10 <- plot_bmd_confidence_intervals(head(bmd_bounds_df, 20))
ggsave(paste0(save_path,"plot10.jpeg"), plot = p10, width = PLOT_WIDTH, height = PLOT_HEIGHT + 500, units = "px", dpi = 600)

