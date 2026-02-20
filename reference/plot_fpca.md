# Plot FPCAs

This function plots the first two principal components calculated

## Usage

``` r
plot_fpca(object, color = NULL)

# S4 method for class 'DGrowthR'
plot_fpca(object, color = NULL)
```

## Arguments

- object:

  DGrowthR object.

- color:

  Covariate by which to color the points in the scatterplot. Should be
  one of the columns in `metadata`. Defaults to NULL

## Value

A scatterplot of the first two FPCs

## Details

\#——————————————————————————————————————
