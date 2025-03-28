% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_metadata}
\alias{example_metadata}
\title{Example DGrowthR Metadata}
\format{
A data frame with 150 rows and 5 variables:
\describe{
\item{curve_id}{A character string indicating the ID or name of the curve.}
\item{Treatment}{A character string indicating the type and, the concentration of the treatment.}
\item{Conc.}{A numeric value indicating the concentration of the treatment.}
\item{Product}{A character string describing the type of product or compound related to the treatment. For example, "Cancer Medication", "Water", "DMSO", or "HIV Medication".}
\item{abx}{A character string indicating whether an antibiotic (or other relevant substance) is present. Values are either "abx" or "no_abx".}
}
}
\usage{
example_metadata
}
\description{
A dataset containing metadata information related to different treatments and their properties.
}
\examples{
head(DGrowthR::example_metadata)

}
\keyword{datasets}
