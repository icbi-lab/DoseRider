% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heatmap.R
\name{dose_response_heatmap}
\alias{dose_response_heatmap}
\title{Create Complex Heatmap for Gene Sets}
\usage{
dose_response_heatmap(
  dose_rider_results,
  dose_col = "Dose",
  dose_unit = "M",
  top = 15,
  order_column = "NegLogPValue",
  decreasing = FALSE,
  fontsize = 6
)
}
\arguments{
\item{dose_rider_results}{A list containing the results of the DoseRider analysis for each gene set.
Each element of the list is a sublist with various metrics and the raw values (expression data) for a gene set.}

\item{dose_col}{A character string specifying the name of the dose column in the raw expression data.}

\item{dose_unit}{A character string specifying the dose units to plot in the column title..}

\item{top}{An integer specifying the number of top gene sets to include in the heatmap. Default is 15.}

\item{order_column}{A character string specifying the column to use for ordering gene sets in the heatmap.}

\item{fontsize}{Integer for fontsize. Default 6}
}
\value{
An object of class `Heatmap` representing the constructed heatmap.
}
\description{
This function creates a heatmap using the `ComplexHeatmap` package, where each row represents a gene set
and each column represents the average expression of genes within that gene set for each dose.
}
\examples{
\dontrun{
  # Assuming `dose_rider_results` is your DoseRider analysis result
  heatmap_plot <- create_complex_heatmap(dose_rider_results, "dose")
  draw(heatmap_plot) # To display the heatmap
}

}
