#' Create a SummarizedExperiment object from assay and metadata data frames
#'
#' This function creates a SummarizedExperiment object which combines the assay
#' data (counts, or other measurements) and associated metadata.
#'
#' @param assay_df Data frame containing assay data, with rows as genes and columns as samples.
#' @param metadata_df Data frame containing metadata for the samples.
#'
#' @return A SummarizedExperiment object with the given assay and metadata.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' assay_df <- data.frame(matrix(runif(100), nrow=10))
#' metadata_df <- data.frame(condition = rep(c("A", "B"), each = 5))
#' se <- create_summarized_experiment(assay_df, metadata_df)
#' @export
create_summarized_experiment <- function(assay_df, metadata_df) {
  # Convert metadata and assay data frames to appropriate objects
  coldata <- as.data.frame(metadata_df)
  assay <- as.matrix(assay_df)

  # Create SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(counts = assay), colData = coldata)

  return(se)
}


#' Prepare data for use with glmmTMB from a SummarizedExperiment or ExpressionSet object
#'
#' This function prepares gene expression data from a `SummarizedExperiment` or `ExpressionSet` object
#' into a long-format data frame suitable for use with the `glmmTMB` function. It includes transformation
#' to long format, integration of metadata, optional log10 transformation of dose values, and generation
#' of spline basis functions for dose-response modeling. When using RNA-seq count data from a
#' `SummarizedExperiment`, it optionally estimates size factors and dispersion parameters.
#'
#' @param se A `SummarizedExperiment` or `ExpressionSet` object containing assay and sample metadata.
#' @param geneset A character vector of gene identifiers to include.
#' @param dose_col Character, name of the dose column in the sample metadata.
#' @param sample_col Character, name of the sample column identifying each observation.
#' @param omic Character string indicating the type of omics data. Set to `"rnaseq"` for count data.
#'             Default is `"rnaseq"`.
#' @param log_transform Logical. Whether to log10-transform dose values. Default is `FALSE`.
#' @param spline_knots Integer. Number of internal spline knots. Must be smaller than the number of unique doses. Default is `2`.
#' @param knot_method Character. Method for selecting spline knot positions. Options:
#'   - `"quantile"` (default): knots placed at quantiles of the dose.
#'   - `"geometric"`: knots spaced in log scale.
#'   - `"manual"`: uses user-specified `manual_knots`.
#' @param manual_knots Numeric vector. Knot positions to use if `knot_method = "manual"`. Ignored otherwise.
#'
#' @return A list with:
#'   - `long_df`: A data frame in long format including expression values, sample metadata, and spline columns.
#'   - `spline_info`: A list with:
#'       - `knots`: Vector of selected internal knot positions.
#'       - `boundary_knots`: The boundary knots used in the spline function.
#'
#' @importFrom reshape melt
#' @importFrom splines bs
#' @examples
#' # Prepare using default quantile-based spline knots
#' data_prep <- prepare_data(se, geneset, dose_col = "dose", sample_col = "sample")
#'
#' # Geometric spacing of spline knots
#' data_prep <- prepare_data(se, geneset, dose_col = "dose", sample_col = "sample", knot_method = "geometric")
#'
#' # Manual knot specification
#' manual_knots <- c(0.01, 0.1, 1, 10)
#' data_prep <- prepare_data(se, geneset, dose_col = "dose", sample_col = "sample", knot_method = "manual", manual_knots = manual_knots)
#'
#' @export
prepare_data <- function(se, geneset, dose_col, sample_col, omic = "rnaseq",
                         log_transform = FALSE, spline_knots = 2, knot_method = "quantile", manual_knots = NULL) {

  # Determine class and extract expression + metadata
  if (inherits(se, "SummarizedExperiment")) {
    expr_data <- assay(se)
    meta_data <- as.data.frame(colData(se))
  } else if (inherits(se, "ExpressionSet")) {
    expr_data <- exprs(se)
    meta_data <- as.data.frame(pData(se))
  } else {
    stop("Unsupported object type. Only SummarizedExperiment or ExpressionSet are allowed.")
  }

  common_genes <- intersect(geneset, rownames(expr_data))

  if (length(common_genes) > 3) {
    long_df <- suppressWarnings(as.data.frame(reshape::melt(t(expr_data[common_genes,]), as.is = TRUE), warning = FALSE))
    colnames(long_df) <- c(sample_col, "gene", "counts")

    if (omic == "rnaseq" && inherits(se, "SummarizedExperiment")) {
      se <- estimate_model_parameters(se)
      long_df$size_factor <- colData(se)$size_factors[match(long_df[[sample_col]], rownames(colData(se)))]
      long_df$theta <- rowData(se)$theta[match(long_df$gene, rownames(rowData(se)))]
      long_df$dispersion <- rowData(se)$dispersions[match(long_df$gene, rownames(rowData(se)))]
    }

    meta_data$sample <- rownames(meta_data)
    meta_data[[sample_col]] <- rownames(meta_data)

    long_df <- merge(long_df, meta_data, by = sample_col)

    long_df[[dose_col]] <- as.numeric(as.character(long_df[[dose_col]]))

    if (log_transform) {
      long_df[[paste0("log_", dose_col)]] <- log10(long_df[[dose_col]] + 1e-6)
      dose_col <- paste0("log_", dose_col)
    }

    unique_doses <- unique(long_df[[dose_col]])
    nk <- min(spline_knots, length(unique_doses) - 1)

    if (nk < 1) {
      knot_positions <- NULL
    } else {
      if (knot_method == "quantile") {
        knot_positions <- stats::quantile(unique_doses, probs = seq(1, nk) / (nk + 1))
      } else if (knot_method == "geometric") {
        if (log_transform) {
          raw_doses <- 10^unique_doses
          raw_knots <- exp(seq(log(min(raw_doses[raw_doses > 0])), log(max(raw_doses)), length.out = nk))
          knot_positions <- log10(raw_knots)
        } else {
          raw_doses <- unique_doses
          knot_positions <- exp(seq(log(min(raw_doses[raw_doses > 0])), log(max(raw_doses)), length.out = nk))
        }
      } else if (knot_method == "manual" && !is.null(manual_knots)) {
        knot_positions <- manual_knots
      } else {
        stop("Invalid knot_method. Choose 'quantile', 'geometric', or 'manual' with 'manual_knots'.")
      }
    }

    spline_dose <- splines::bs(long_df[[dose_col]], knots = knot_positions,
                               degree = 3,
                               Boundary.knots = range(long_df[[dose_col]]), intercept = FALSE)

    spline_df <- as.data.frame(spline_dose)
    colnames(spline_df) <- paste0("CubicSpline_", seq_len(ncol(spline_df)))

    long_df <- cbind(long_df, spline_df)
    long_df <- as.data.frame(long_df)
  } else {
    long_df <- NULL
    spline_dose <- NULL
  }

  return(list(long_df = long_df, spline_info = spline_dose))
}





#' Update SummarizedExperiment with Estimated Model Parameters
#'
#' Estimates necessary parameters for DESeq2 or edgeR and updates the SummarizedExperiment object
#' with size factors, dispersions, and theta values in its rowData and colData.
#'
#' @param se SummarizedExperiment object containing count data and sample metadata.
#' @return Updated SummarizedExperiment object with size factors, dispersions, and theta values included.
#' @importFrom edgeR DGEList calcNormFactors estimateDisp
#' @export
estimate_model_parameters <- function(se) {

    # Create DGEList object
    dge <- DGEList(counts = assay(se), group = se$dose)

    # Estimate normalization factors and dispersion
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge)

    # Add size factors to colData
    sizeFactors <- dge$samples$norm.factors * dge$samples$lib.size
    se$size_factors <- sizeFactors

    # Add dispersions to rowData
    dispersions <- dge$common.dispersion
    theta <- 1 / dispersions
    rowData(se)$dispersion <- dispersions
    rowData(se)$theta <- theta
    colData(se)$size_factors <- sizeFactors


  return(se)
}
