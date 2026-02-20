# Plot UMAP components

This function plots the first two UMAP components calculated

## Usage

``` r
plot_umap(object, color = NULL)

# S4 method for class 'DGrowthR'
plot_umap(object, color = NULL)
```

## Arguments

- object:

  Degrowth object.

- color:

  Covariate by which to color the points in the scatterplot. Should be
  one of the columns in `metadata`. Defaults to NULL

## Value

A scatterplot of the first two UMAP components
