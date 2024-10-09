# DoseRider: A multi-omics approach to study dose-response relationships at the pathway level using mixed models

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/doseRider)](https://cran.r-project.org/package=doseRider)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Overview

DoseRider enhances toxicogenomics by employing mixed models with cubic splines to study dose-response relationships at the pathway level. This approach is ideal for multi-omics research, supporting both RNA-seq and metabolomics data, and is available as an R package and [web application](https://doserider.i-med.ac.at/). DoseRider identifies pathway trends, calculates Trend Change Doses (TCD), and estimates Benchmark Doses (BMD). This enables a deeper understanding of the molecular responses to compounds or drugs across various concentrations.

## Key Features

- **Non-linear and Linear Modeling:** Supports linear mixed models with cubic splines and generalized mixed models to handle complex dose-response relationships.
- **Pathway-Level Dose-Response Modeling:** Models dose-response at the pathway or gene-set level, calculating TCDs and BMDs.
- **Multi-Omics:** Compatible with various omics data, including RNA-Seq and metabolomics.
- **Parallel Computing:** Utilizes parallel processing to efficiently handle large datasets.
- **Visualization Tools:** Multiple built-in plotting functions to visualize dose-response trends and model outputs.
- **Customizable Gene Sets:** Filter and analyze custom or preprocessed gene sets.

## Installation

To install the latest development version from GitHub, use the following:

```r
# install.packages("devtools")
devtools::install_github("icbi-lab/doseRider")
```

## Usage Example

Here is an example of using DoseRider to analyze RNA-Seq data:

```r
# Load DoseRider
library(doseRider)

# Load your gene expression data
data("bpaf_data")

# Load gene sets
gmt_path <- system.file("extdata", "High-Response-Toxicogenomics.gmt", package = "doseRider")
gmt <- read_gmt(gmt_path)

# Perform dose-response analysis
# Run doseRider analysis
dose_rider_results <- DoseRiderParallel(
  se = bpaf_data, 
  gmt = gmt, 
  dose_col = "Dose", 
  omic = "rnaseq", 
  minGSsize = 20, 
  maxGSsize = 200, 
  method = "fdr", 
  covariates = c(),
  modelType = "LMM", 
  num_cores = 10,
  FilterPathway = TRUE,
  log_transform = TRUE,
  models = c("linear", "non_linear_mixed")
)
```

### Visualization

DoseRider provides various functions to visualize dose-response data. Below are some key visualizations:

1. **Dose-Response Heatmap:**

```r
p1 <- dose_response_heatmap(dose_rider_results_filter, dose_col = "Dose", top = 15)
jpeg(file=paste0(save_path, "plot1.jpeg"), width = 10, height = 10, units = "cm", res = 600)
plot(p1, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
```

2. **Gene Set Random Effects Plot:**

```r
p2 <- plot_gene_set_random_effects(dose_rider_results_filter, dose_col = "log_Dose", top = 15)
ggsave(paste0(save_path, "plot2.jpeg"), plot = p2, width = 10, height = 10, units = "cm", dpi = 600)
```

3. **Top Pathway Responses:**

```r
p3 <- plot_top_pathway_responses(dose_rider_results_filter, top = 2, ncol = 2, text_size = 5, dose_col = "log_Dose", clusterResults = TRUE)
ggsave(paste0(save_path, "plot3.jpeg"), plot = p3, width = 10, height = 10, units = "cm", dpi = 600)
```

4. **Gene Random Effect Relationship Plot:**

```r
p4 <- plot_gene_random_effect_relationship(dose_rider_results_filter, "Estrogen signaling pathway - Homo sapiens (human)")
ggsave(paste0(save_path, "plot4.jpeg"), plot = p4, width = 10, height = 10, units = "cm", dpi = 600)
```

5. **Dot Plot of Top Pathways:**

```r
p5 <- plot_dotplot_top_pathways(dose_rider_results_filter, top = 15)
ggsave(paste0(save_path, "plot5.jpeg"), plot = p5, width = 10, height = 10, units = "cm", dpi = 600)
```

6. **Gene Heatmap for a Specific Pathway:**

```r
p6 <- create_gene_heatmap(dose_rider_results_filter, dose_col = "Dose", gene_set_name = "Estrogen signaling pathway - Homo sapiens (human)")
jpeg(file=paste0(save_path,"plot6.jpeg"), width = 10, height = 10, units = "cm", res = 600)
plot(p6, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
```

### BMD and TCD Calculations

DoseRider allows the calculation of BMD and TCD values:

```r
# Calculate BMD range and plot the density and peaks
data_bmd <- get_bmd_range(dose_rider_results = dose_rider_results_filter)
p7 <- plot_bmd_density_and_peaks(data_bmd)
ggsave(paste0(save_path, "plot7.jpeg"), plot = p7, width = 10, height = 10, units = "cm", dpi = 600)

# Plot BMD confidence intervals
p10 <- plot_bmd_confidence_intervals(head(bmd_bounds_df, 20))
ggsave(paste0(save_path, "plot10.jpeg"), plot = p10, width = 10, height = 10, units = "cm", dpi = 600)
```

## Toxicogenomics Gene Set

The **Toxicogenomics Gene Set** focuses on pathways with significant changes across compounds in the TG-GATES database. The score is calculated by multiplying NES with the negative logarithm of the p-value per dose level, and then averaged across doses.

![Gene Heatmap for Specific Pathway](./plots/HeatmapHighResponsive.jpeg)

## Contributing

We welcome contributions! Open an issue or submit a pull request on our [GitHub repository](https://github.com/icbi-lab/doserider).

## License

DoseRider is licensed under the MIT License. For more information, see the [LICENSE](LICENSE) file.
