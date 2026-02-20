# Generate predictions from a Gaussian Process model

This function generates predictions from a Gaussian Process model fitted
to an object. The function supports "alternative" and "null" models. The
"alternative" model predicts based on both time and contrast, while the
"null" model predicts based on time only.

## Usage

``` r
.predict_gp(
  XX,
  lagp_fit,
  complete_sigma = FALSE,
  prepare_dataframe = FALSE,
  estimate_derivatives = TRUE,
  delete = TRUE
)
```

## Arguments

- XX:

  The covariate matrix for which we want ot make predictions

- lagp_fit:

  The output of an laGP fit.

- complete_sigma:

  A logical value indicating whether the full covariance matrix should
  be returned. This should be set to TRUE if sampling the posterior is
  desired.

- estimate_derivatives:

  A logical value indicating if the empirical first and second
  derivatives should be returned.

- delete:

  A logical value. If TRUE, the Gaussian Process model will be deleted
  from the object after prediction.

## Value

The mean predicted values.

## See also

[`laGP::predGPsep()`](https://rdrr.io/pkg/laGP/man/predGP.html)
