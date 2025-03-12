# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased] - 2025-03-12  

### **Added**  
- **Enhanced Knot Selection in `prepare_data()` for `DoseRider()`**  
  - Introduced support for selecting the number of knots (`nk`).  
  - Implemented three **knot selection methods**:  
    - **Quantile-Based (`"quantile"`)**: Places knots at evenly spaced quantiles of unique dose values.  
    - **Geometric Progression (`"geometric"`)**: Places knots based on an exponential scale for log-transformed doses.  
    - **Manual (`"manual"`)**: Allows users to specify custom knot positions.  

### **Fixed**  
- Resolved minor bugs affecting knot selection in `prepare_data()`.  
- Improved stability when handling log-transformed doses with geometric knot placement.  
- Fixed inconsistencies in function documentation.



## [Unreleased] - 2025-03-05

### Added
- **New Modular Pathway Filtering System**
  - Introduced `apply_pathway_filters()` to dynamically apply selected filters based on user-defined parameters.
  - **PCA-Based Filtering (`"pca"`)**: Identifies and removes pathways with strong antagonist expression patterns.
  - **Variance-Based Filtering (`"variance"`)**: Removes gene sets with low expression variability.
  - **Correlation-Based Filtering (`"correlation"`)**: Filters pathways with a high proportion of negatively correlated genes.
  - **Customizable Filtering Parameters**: Users can now define thresholds for PCA variance (`pca_threshold`), variance percentile (`variance_percentile`), and correlation cutoff (`correlation_threshold`).

- **New Parameter: `FilterPathway` in `DoseRider()`**
  - Enables or disables pathway filtering.
  - Allows selection of multiple filtering methods.

- **Performance Improvements**
  - Optimized correlation computation by pivoting data before pairwise correlation analysis.
  - Improved handling of log-transformed dose values in model fitting.
  
### Fixed
- Resolved issue where PCA filtering incorrectly excluded pathways with moderate variance.
- Fixed incorrect data pivoting before computing correlation matrices.
- Corrected function documentation inconsistencies.

---

## [Unreleased] - 2024-10-27

### Added
- **New functions:** `doseRiderLMM` and `doseRiderGAMM`
  - **Description:** Introduced support for Linear Mixed Models (LMM) and Generalized Additive Mixed Models (GAMM).

- **Enhancements in `lmm.R`**
  - **New Functions:** `fit_lmm` and `create_lmm_formula`
    - `fit_lmm`: Fits linear mixed models.
    - `create_lmm_formula`: Constructs LMM formulas dynamically.
  - **New Parameter:** `omic`
    - Allows users to choose between Gaussian and negative binomial models.

- **New function: `smooth_pathway_trend`**
  - **Description:** Uses `predict()` with `re.form = NA` to smooth pathway trends.

- **Initial Creation of `bmd.R`**
  - **Purpose:** Benchmark Dose (BMD) estimation (still in development).

---

## [Unreleased] - 2024-08-03

### Added
- **New option:** `center_values` in `plot_smooth`
  - **Description:** Enables centering and scaling expression values per gene.

- **New function:** `DoseRiderParallel`
  - **Description:** Allows parallel processing of gene sets using multiple cores for improved computation speed.
  - **Usage:**  
    ```r
    DoseRiderParallel(se, gmt, dose_col = "dose", sample_col = "sample", omic = "rnaseq", minGSsize = 5, maxGSsize = 300, method = "fdr", num_cores = 5)
    ```

### Changed
- **Improved Calculation of `p_value_cubic` in `DoseRider()`**
  - Now computes `p.value` for both linear and cubic splines.

---

## [Unreleased] - 2024-05-01

### Added
- **Initial Release of DoseRider**
- **Implemented GLMM with Splines** for dose-response analysis.
- **Web-based Visualization Interface**
- **Support for RNA-seq and Other Omics Data**

---


---

For bug reports and contributions, please visit our [GitHub repository](https://github.com/icbi-lab/DoseRider).  

