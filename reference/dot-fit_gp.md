# Fit a Gaussian process model to data

This function fits a Gaussian process model to the provided data.

## Usage

``` r
.fit_gp(X_data, delete = FALSE)
```

## Arguments

- X_data:

  A data.frame containing the covariates used for model fitting. Should
  have at least timepoint, timepoint_n, and od columns.

- delete:

  A logical value. If TRUE, the Gaussian process model will be deleted
  after fitting.

## Value

The input object, with the Gaussian process model information as a list
