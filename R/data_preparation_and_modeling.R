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


#' Prepare data for use with glmmTMB from a SummarizedExperiment object
#'
#' This function prepares the data from a SummarizedExperiment object into a format
#' suitable for the `glmmTMB` function. It transforms the data into long format,
#' adds size factors and dispersions (if applicable), applies log10 transformation to doses,
#' and creates splines for dose using different knot selection methods.
#'
#' @param se A `SummarizedExperiment` object containing both assay and metadata.
#' @param geneset A character vector of gene names to be considered.
#' @param dose_col The name of the dose column in the metadata of the `SummarizedExperiment` object.
#' @param sample_col The name of the sample column in the metadata of the `SummarizedExperiment` object.
#' @param omic A character string specifying the type of omics data. Default is `"rnaseq"`.
#' @param log_transform Logical, whether to log10 transform the dose values. Default is `FALSE`.
#' @param spline_knots Integer, the number of internal knots for the spline function. Default is `2`.
#'                     The number of knots should be lower than the number of unique dose values.
#' @param knot_method Character, the method for selecting spline knots. Options are:
#'   - `"quantile"` (default): Selects knots based on percentiles of the dose distribution.
#'   - `"geometric"`: Places knots at logarithmically spaced intervals.
#'   - `"manual"`: Uses user-defined knots from the `manual_knots` parameter.
#' @param manual_knots Numeric vector, custom knot positions (used only if `knot_method = "manual"`).
#'                     If `NULL`, this option is ignored.
#'
#' @return A list containing:
#'   - `long_df`: A data frame in long format with precomputed spline columns (`CubicSpline_1`, etc.).
#'   - `spline_info`: A list with:
#'       - `knots`: The selected knot positions.
#'       - `boundary_knots`: The boundary knot positions.
#'
#' @importFrom reshape melt
#' @importFrom splines bs
#' @examples
#' # Example using quantile-based knots (default)
#' data_prep <- prepare_data(se, geneset, dose_col, sample_col)
#'
#' # Example using geometric knots (log-spaced)
#' data_prep <- prepare_data(se, geneset, dose_col, sample_col, knot_method = "geometric")
#'
#' # Example using manually selected knots
#' manual_knot_positions <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 10)
#' data_prep <- prepare_data(se, geneset, dose_col, sample_col, knot_method = "manual", manual_knots = manual_knot_positions)
#' @export
prepare_data <- function(se, geneset, dose_col, sample_col, omic = "rnaseq",
                         log_transform = FALSE, spline_knots = 2, knot_method = "quantile", manual_knots = NULL) {

    # Prepare data for glmmTMB, including size factors and dispersions
    common_genes <- intersect(geneset, rownames(se))

    if (length(common_genes) > 3) {
      long_df <- suppressWarnings(as.data.frame(reshape::melt(t(assay(se)[common_genes,]), as.is = TRUE), warning = FALSE))
      colnames(long_df) <- c(sample_col, "gene", "counts")

      if (omic == "rnaseq") {
        se <- estimate_model_parameters(se)
        long_df$size_factor <- colData(se)$size_factors[match(long_df[[sample_col]], rownames(colData(se)))]
        long_df$theta <- rowData(se)$theta[match(long_df$gene, rownames(rowData(se)))]
        long_df$dispersion <- rowData(se)$dispersions[match(long_df$gene, rownames(rowData(se)))]
      }

      colData(se)$sample <- colnames(se)
      colData(se)[[sample_col]] <- colnames(se)

      long_df <- merge(long_df, colData(se), by = sample_col)

      # Convert dose to numeric
      long_df[[dose_col]] <- as.numeric(as.character(long_df[[dose_col]]))

      # Log10 transform doses if specified
      if (log_transform) {
        long_df[[paste0("log_", dose_col)]] <- log10(long_df[[dose_col]] + 1e-6)  # Avoid log(0)
        dose_col <- paste0("log_", dose_col)
      }

      # --- CONTROLLED SPLINE KNOT SELECTION ---
      unique_doses <- unique(long_df[[dose_col]])
      nk <- min(spline_knots, length(unique_doses) - 1)  # Ensure no excess knots

      if (nk < 1) {
        knot_positions <- NULL  # No internal knots if not enough unique doses
      } else {
        if (knot_method == "quantile") {
          knot_positions <- stats::quantile(unique_doses, probs = seq(1, nk) / (nk + 1))
        } else if (knot_method == "geometric") {
          if (log_transform){
            # Undo log10 transformation: Convert log-doses back to original scale
            raw_doses <- 10^unique_doses

            # Select knots in raw (non-log) scale using geometric progression
            raw_knots <- exp(seq(log(min(raw_doses[raw_doses > 0])),
                                 log(max(raw_doses)), length.out = nk))

            # Convert knots back to log10 scale
            knot_positions <- log10(raw_knots)
          } else {
            raw_doses <- unique_doses
            knot_positions <- exp(seq(log(min(raw_doses[unique_doses > 0])),
                                      log(max(unique_doses)), length.out = nk))

          }
        } else if (knot_method == "manual" && !is.null(manual_knots)) {
          knot_positions <- manual_knots
        } else {
          knot_positions <- 0
          stop("Invalid knot_method. Choose 'quantile', 'geometric', or 'manual' with 'manual_knots'.")
        }
      }

      # Generate B-splines
      spline_dose <- splines::bs(long_df[[dose_col]], knots = knot_positions,
                                 degree = 3,
                                 Boundary.knots = range(long_df[[dose_col]]), intercept = FALSE)

      # Convert to DataFrame
      spline_df <- as.data.frame(spline_dose)
      colnames(spline_df) <- paste0("CubicSpline_", seq_len(ncol(spline_df)))

      # Merge with long_df
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
