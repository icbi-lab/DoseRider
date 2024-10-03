#' Compute IC50 from SummarizedExperiment object
#'
#' This function computes the IC50 (half-maximal inhibitory concentration) for each gene
#' in a given `SummarizedExperiment` object based on the dose-response relationship.
#' The function takes the name of the dose column and the viability column (response)
#' from the `colData` and returns a vector of IC50 values.
#'
#' @param se_object A `SummarizedExperiment` object containing the dose-response data.
#' @param dose_col A string representing the name of the column in `colData` that contains the dose information.
#' @param viability_col A string representing the name of the column in `colData` that contains the viability or response data.
#'
#' @return A named vector of IC50 values for each gene, where the names are the gene identifiers.
#'
#' @examples
#' # Assuming `se` is a SummarizedExperiment object with dose and viability information
#' # ic50_values <- compute_IC50(se, dose_col = "dose", viability_col = "viability")
#'
#' @import SummarizedExperiment
#' @import drc
#' @export
compute_IC50 <- function(se_object, dose_col, viability_col) {

  # Extract expression matrix and metadata from the SummarizedExperiment object
  expression_data <- assay(se_object)
  metadata <- colData(se_object)

  # Check if dose and viability columns are present in metadata
  if (!(dose_col %in% colnames(metadata))) {
    stop(paste("The metadata must contain a", dose_col, "column."))
  }

  if (!(viability_col %in% colnames(metadata))) {
    stop(paste("The metadata must contain a", viability_col, "column."))
  }

  # Initialize a vector to store IC50 values for each gene
  ic50_results <- numeric(length(rownames(expression_data)))
  names(ic50_results) <- rownames(expression_data)

  # Iterate over each gene (row) in the expression data
  for (i in seq_len(nrow(expression_data))) {

    # Get expression values for the current gene across samples
    gene_expression <- expression_data[i, ]

    # Create a data frame with dose and viability (response) values
    data <- data.frame(
      dose = as.numeric(metadata[[dose_col]]),
      response = as.numeric(metadata[[viability_col]])
    )

    # Fit a dose-response curve using the drc package (log-logistic model)
    fit <- tryCatch({
      drm(response ~ dose, data = data, fct = LL.4())
    }, error = function(e) {
      NA  # Handle cases where fitting fails
    })

    # Extract the IC50 value if fitting was successful
    if (!is.na(fit)) {
      ic50_results[i] <- ED(fit, 50)["estimate"]
    } else {
      ic50_results[i] <- NA  # Set to NA if fitting failed
    }
  }

  # Return the vector of IC50 values
  return(ic50_results)
}

