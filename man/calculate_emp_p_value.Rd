% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculations.R
\name{calculate_emp_p_value}
\alias{calculate_emp_p_value}
\alias{calculate_emp_p_value,DGrowthR-method}
\title{Calculate the Empirical p_value}
\usage{
calculate_emp_p_value(object)

\S4method{calculate_emp_p_value}{DGrowthR}(object)
}
\arguments{
\item{object}{A DGrowthR object which contains the Gaussian Process regression results.}
}
\value{
A numeric value representing the calculated Empirical p-value.
}
\description{
This function calculates the Empirical p-value (emp_p_value) based on a given object.
}
