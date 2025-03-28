% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/growth_parameters.R
\name{.characterize_growth}
\alias{.characterize_growth}
\title{Estimate growth parameters}
\usage{
.characterize_growth(df, od_auc_at_t = NULL)
}
\arguments{
\item{df}{A data.frame with timepoint and mean fields.}

\item{od_auc_at_t}{A numeric value for which od and AUC are calculated for}
}
\value{
A data.frame with the requested growth parameters for the curve given as input
}
\description{
Estimate growth parameters
}
\keyword{internal}
