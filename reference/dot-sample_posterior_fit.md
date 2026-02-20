# Sample the posterior GP model for new curves

Sample the posterior GP model for new curves

## Usage

``` r
.sample_posterior_fit(lagp_fit, n_samples = 100)
```

## Arguments

- lagp_fit:

  The output of an laGP fit.

- n_samples:

  A numberic value indicating the number of curves to sample from
  posterior

## Value

A data.frame with sample_id, timepoint, and mean columns.

## See also

[`laGP::predGPsep()`](https://rdrr.io/pkg/laGP/man/predGP.html)
