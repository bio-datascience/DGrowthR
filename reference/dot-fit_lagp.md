# Fit a Gaussian process model to data with laGP

This function fits a Gaussian process model to the provided data. It is
an internal function that really only does the fit. Nothing else is done
or returned by this function.

## Usage

``` r
.fit_lagp(X, y, delete = FALSE)
```

## Arguments

- X:

  A data.frame containing the covariates used for model fitting.

- y:

  A vector containing the values to be predicted.

- delete:

  A logical value. If TRUE, the Gaussian process model will be deleted
  after fitting.

## Value

The input object, with the Gaussian process information as a list.
