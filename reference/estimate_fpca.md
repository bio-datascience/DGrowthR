# Functional Principal Component Analysis

This function estimates the functional principal components of the
complete dataset using the fdapace package

## Usage

``` r
estimate_fpca(object, data = "od_data", ...)

# S4 method for class 'DGrowthR'
estimate_fpca(object, data = "od_data", ...)
```

## Arguments

- object:

  DGrowthR object.

- data:

  A string indicating the source data to embed. At the moment only
  embeds "od_data"

- ...:

  Additional arguments to pass to the optns argument of FPCA

## Value

The result from FPCA function from fdapace package

## See also

[`fdapace::FPCA()`](https://rdrr.io/pkg/fdapace/man/FPCA.html)
