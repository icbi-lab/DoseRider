% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmd.R
\name{get_bmd_range}
\alias{get_bmd_range}
\title{Compute BMD Range and Find Peaks}
\usage{
get_bmd_range(dose_rider_results)
}
\arguments{
\item{dose_rider_results}{A list containing the results of DoseRider analysis.}
}
\value{
A list containing the BMD values (`x`), their density (`y`), and the identified peaks (`bmd`).
}
\description{
This function computes the range of Benchmark Dose (BMD) values from the DoseRider results and identifies peaks in the density of BMD values.
}
\examples{
\dontrun{
dose_rider_results <- DoseRider(se, gmt, "dose", "sample")
bmd_range <- get_bmd_range(dose_rider_results)
print(bmd_range)
}

}
