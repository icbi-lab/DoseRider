# doseRider

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/doseRider)](https://cran.r-project.org/package=doseRider)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Overview

doseRider is an advanced R package designed for comprehensive analysis of dose-response relationships in gene expression data. It employs both Linear Mixed Models with cubic splines and Generalized Mixed Models, making it adept at handling complex, non-linear dose-response patterns. This package is tailored for a multi-omic approach, seamlessly adapting to various types of omic data, including RNA-Seq and microarray.

## Key Features

- **Versatile Modeling Capability:** Implements Linear Mixed Models with cubic splines and Generalized Mixed Models to accommodate non-linear dose-response relationships.
- **Breakpoint Identification:** Extract significant breakpoints in dose-response curves to identify critical dose levels impacting gene expression.
- **Multi-Omic Approach:** Automatically adjusts to different omic technologies, choosing appropriate statistical models for RNA-Seq or microarray data.
- **Automatic Model Selection:** Dynamically selects the most significant model between null, linear, and non-linear models, ensuring robust analysis.
- **Omics-Specific Distribution Selection:** Automatically opts for Gaussian or negative binomial distribution based on the omic technology used, enhancing the accuracy of the analysis.
- **Parallel Processing:** Utilizes parallel computing for efficient processing of large-scale datasets.

## Installation

Or install the development version from GitHub:

```R
# install.packages("devtools")
devtools::install_github("icbi-lab/doseRider")
```

## Usage Example

Below is an example workflow demonstrating the use of doseRider:

```R
# Load the doseRider package
library(doseRider)

# Example: Reading gene expression data as a SummarizedExperiment object
data <- readRDS("gene_expression_data.rds")

# Load gene set collections
gmt <- loadCPDB("Symbol")

# Filter gene sets by size
gmt_filtered <- filter_gmt_by_size(gmt, minGenesetSize = 10, maxGenesetSize = 50)

# Perform dose-response analysis
result <- DoseRiderParallel(se = data, gmt = gmt_filtered, dose_col = "dose", 
                            sample_col = "sample", omic = "rnaseq", minGSsize = 10, 
                            maxGSsize = 200, num_cores = 4)

```

### Visualization

Visualizations in doseRider provide insightful representations of the dose-response data. The package offers various plotting functions to explore and interpret the results effectively. Here's an overview of some key visualizations:

1. **Dose Response Heatmap:**

```{r}
   p1 <- dose_response_heatmap(dose_rider_results, dose_col = "Dose", top = 10)
 ```
 This heatmap displays the average gene expression response across different doses for the top 10 gene sets. Each row represents a gene set, and each column corresponds to a specific dose. The color intensity in each cell reflects the magnitude of gene expression, providing a clear visual representation of how gene sets respond to varying doses.

2. **Gene Set Random Effects Plot:**

```{r}
   p2 <- plot_gene_set_random_effects(dose_rider_results, dose_col = "Dose", top = 10)
```
   
The plot showcases the distribution of random effects for genes within each of the top 10 gene sets. This visualization helps in understanding the variability of gene expression within each gene set. A wider distribution indicates greater variability, offering insights into the heterogeneity of gene responses within the set.

3. **Top Pathway Response Plots:**

```{r}
   p3 <- plot_top_pathway_responses(dose_rider_results, top = 10, ncol = 5)
```
  
  This function generates dose-response plots for the top 10 pathways, arranged in a grid layout. Each plot illustrates the dose-response relationship for a specific gene set, highlighting the trends and patterns of gene expression as a function of the dose. These plots are instrumental in identifying gene sets with significant dose-dependent changes.

## Contributing

We welcome contributions to doseRider! If you find bugs, have feature requests, or want to contribute to the development, please open an issue or submit a pull request on our [GitHub repository](https://github.com/icbi-lab/doseRider).

## License

doseRider is licensed under the MIT License. For more details, see the [LICENSE](LICENSE) file.

