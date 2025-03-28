% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_likelihood_ratio}
\alias{plot_likelihood_ratio}
\alias{plot_likelihood_ratio,DGrowthR-method}
\title{Plot Likelihood Ratios}
\usage{
plot_likelihood_ratio(object, num_bins = 50)

\S4method{plot_likelihood_ratio}{DGrowthR}(object, num_bins = 50)
}
\arguments{
\item{object}{An object of class "DGrowthR".}

\item{num_bins}{Number of bins to be used in the histogram.}
}
\value{
A ggplot object representing the histogram plot.
}
\description{
This function plots a histogram of permuted Likelihood Ratios with the actual Likelihood Ratio indicated by a vertical line.
}
