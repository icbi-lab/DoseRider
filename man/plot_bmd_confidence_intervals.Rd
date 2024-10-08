% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pathway_response.R
\name{plot_bmd_confidence_intervals}
\alias{plot_bmd_confidence_intervals}
\title{Plot Bubble Plot of BMD with Confidence Intervals}
\usage{
plot_bmd_confidence_intervals(bmd_bounds_df, top = 10)
}
\arguments{
\item{bmd_bounds_df}{A data frame containing the BMD results with columns for lower bound, upper bound, mean BMD, and median BMD.}
}
\value{
A ggplot object representing the bubble plot of BMD with confidence intervals.
}
\description{
This function creates a bubble plot visualizing the BMD results with confidence intervals.
The size of the bubbles represents the median BMD, while the error bars show the confidence intervals.
}
\examples{
\dontrun{
  # Assuming bmd_bounds_df is available
  bmd_plot <- plot_bmd_confidence_intervals(bmd_bounds_df, top = 10)
  print(bmd_plot)
}

}
