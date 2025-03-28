% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_gp}
\alias{plot_gp}
\alias{plot_gp,DGrowthR-method}
\title{Plot Gaussian Process}
\usage{
plot_gp(object, gpfit_id = NULL, show_input_od = FALSE)

\S4method{plot_gp}{DGrowthR}(object, gpfit_id = NULL, show_input_od = FALSE)
}
\arguments{
\item{object}{An object of class "DGrowthR"}

\item{gpfit_id}{The name of the GP model that should be plotted. Should be in the gpfit_info$gpfit_data}

\item{show_input_od}{A logical value indicating if the od curves used to create the GP should be plotted}
}
\value{
A ggplot object showing lineplots of growth curves and GP mean.
}
\description{
This function plots a specific gaussian process with the 95\% confidence interval.
}
