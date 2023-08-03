# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased] - 2023-08-03

### Added
- New option: `center_values` in `plot_smooth`.
  - Description: This option allows for centering and scaling the expression values for each gene.
  - Usage: When `center_values` is enabled `T`, the expression values of each gene are centered around their mean and scaled by dividing by their standard deviation.

### Changed
- Calculation of `p_value_cubic` in `DoseRider` for the comparison between linear and cubic splines.
  - Description: Previously, only the linear base vs cubic `p.value` was calculated during the comparison. Now, we also calculate the `p.value` for the cubic spline.
  - Reason: To have obtain the most non-linear significant genesets.

### Fixed
 
