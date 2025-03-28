% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DGrowthR.R
\name{gather_od_data}
\alias{gather_od_data}
\alias{gather_od_data,DGrowthR-method}
\title{Gather specific growth data from the DGrowthR object}
\usage{
gather_od_data(
  object,
  metadata_field,
  field_value,
  downsample_every_n_timepoints = 1
)

\S4method{gather_od_data}{DGrowthR}(
  object,
  metadata_field,
  field_value,
  downsample_every_n_timepoints = 1
)
}
\arguments{
\item{object}{An object of class "DGrowthR".}

\item{metadata_field}{A character string specifying one of the fields in the metadata from which the variables will be searched in.}

\item{field_value}{The value of the metadata_field that for which the relevant growth curves will be extracted.}

\item{downsample_every_n_timepoints}{A numeric value indicating that the OD from every n timepoint should be used for GP fit. Might seep up fitting.}
}
\value{
A data.frame containing the optical density data requested.
}
\description{
Extracts the requested growth curves according to the requested metadata fields and variables.
}
