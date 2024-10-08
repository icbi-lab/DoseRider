% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneset.R
\name{read_gmt}
\alias{read_gmt}
\title{Read GMT File from MSigDB}
\usage{
read_gmt(file_path)
}
\arguments{
\item{file_path}{The path to the GMT file to be read.}
}
\value{
A list where each element is a list with two elements: `$pathway`, the name of the pathway,
and `$genes`, a vector of gene symbols associated with that pathway.
}
\description{
Parses a GMT (Gene Matrix Transposed) file from the Molecular Signatures Database (MSigDB),
structuring it into a list of pathways and their associated genes. Each element of the list
is a list itself, containing the pathway name and a vector of gene symbols.
}
\examples{
file_path <- "path/to/your/c2.cgp.v2023.2.Hs.symbols.gmt"
gmt_data <- read_gmt(file_path)
if (length(gmt_data) > 0) {
  print(gmt_data[[1]])
  print(gmt_data[[2]])
}
}
