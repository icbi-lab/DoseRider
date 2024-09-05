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
#' suitable for the glmmTMB function. This includes transforming the data to long format,
#' adding size factors and dispersions (if applicable), applying log10 transformation to doses,
#' creating splines for dose, and merging with sample metadata.
#'
#' @param se A SummarizedExperiment object containing both assay and metadata.
#' @param geneset A character vector of gene names to be considered.
#' @param dose_col The name of the dose column in the metadata of the SummarizedExperiment object.
#' @param sample_col The name of the sample column in the metadata of the SummarizedExperiment object.
#' @param omic A character string specifying the type of omics data. Default is "rnaseq".
#' @param log_transform Logical, whether to log10 transform the dose values. Default is FALSE.
#' @param spline_knots Number of internal knots to use for splines on dose. Default is 3.
#'
#' @return A data frame in a format suitable for glmmTMB function.
#' @importFrom reshape melt
#' @importFrom splines ns
#'
#' @examples
#' geneset <- rownames(assay_df)
#' dose_col <- "condition"
#' sample_col <- "sample"
#' long_df <- prepare_data(se, geneset, dose_col, sample_col)
#' @export
prepare_data <- function(se, geneset, dose_col, sample_col, omic = "rnaseq", log_transform = FALSE, spline_knots = 3) {

  # Prepare data for glmmTMB, including size factors and dispersions
  common_genes <- intersect(geneset, rownames(se))

  if (length(common_genes) > 0) {
    long_df <- suppressWarnings(as.data.frame(reshape::melt(t(assay(se)[common_genes,]), as.is = TRUE), warning = FALSE))
    colnames(long_df) <- c(sample_col, "gene", "counts")

    if (omic == "rnaseq") {
      # Match size_factors, theta, and dispersions
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
      # Remove rows where dose is 0
      long_df <- long_df[long_df[[dose_col]] > 0, ]
      long_df[[paste0("log_", dose_col)]] <- log10(long_df[[dose_col]])
    }

    # Add natural splines for dose
    dose_splines <- bs(long_df[[dose_col]], knots = spline_knots,
                                              Boundary.knots = range(long_df[[dose_col]]), intercept = FALSE)

    # Store the spline attributes (knots, boundary knots) for later use
    spline_info <- list(knots = attr(dose_splines, "knots"), Boundary.knots = attr(dose_splines, "Boundary.knots"))


    dose_splines_df <- as.data.frame(dose_splines)
    # Multiply spline values by 10, as in the other tool
    dose_splines_df <- dose_splines_df * 10

    # Rename spline columns
    colnames(dose_splines) <- paste("spline_dose_", seq_len(ncol(dose_splines)), sep="")

    # Add splines to the data frame
    long_df <- cbind(long_df, dose_splines)

    long_df <- as.data.frame(long_df, warning = FALSE)
  } else {
    long_df <- NULL
    spline_info <- NULL
  }

  return(list(long_df = long_df, spline_info = spline_info))
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
    sizeFactors <- dge$samples$norm.factors
    se$size_factors <- sizeFactors

    # Add dispersions to rowData
    dispersions <- dge$common.dispersion
    theta <- 1 / dispersions
    rowData(se)$dispersion <- dispersions
    rowData(se)$theta <- theta

  return(se)
}
