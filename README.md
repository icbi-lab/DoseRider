

# doseRider

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/doseRider)](https://cran.r-project.org/package=doseRider)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Overview

The `doseRider` package provides tools for analyzing dose-response relationships in gene expression data using Generalized Additive Mixed Models (GAMMs). It allows users to model and analyze the relationship between gene expression and dose levels of different compounds, identify dose-response patterns, and extract breakpoints. The package is particularly useful for studying the impact of chemical compounds on gene expression profiles and assessing potential toxicological risks.

## Installation

You can install the latest version of `doseRider` from CRAN:

```R
install.packages("doseRider")
```

Or you can install the development version from GitHub:

```R
# install.packages("devtools")
devtools::install_github("username/doseRider")
```

## Usage

```R
# Load the doseRider package
library(doseRider)

# Read in your gene expression data as a SummarizedExperiment object
data <- readRDS("gene_expression_data.rds")

# Perform the dose-response analysis
result <- doseRider(data)

# Extract breakpoints from the results
breakpoints <- extract_breakpoints(result)

# Plot the dose-response curves
plot_dose_response(result)

# Perform downstream analysis and interpretation
# ...
```

## Contributing

Contributions to `doseRider` are welcome. If you find any bugs, have feature requests, or want to contribute improvements or new features, please open an issue or submit a pull request on GitHub.

## License

This package is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
