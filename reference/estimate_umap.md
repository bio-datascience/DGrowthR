# Functional UMAP

This function estimates the UMAP coordinates using the uwot package

## Usage

``` r
estimate_umap(object, data = "od_data", ...)

# S4 method for class 'DGrowthR'
estimate_umap(object, data = "od_data", ...)
```

## Arguments

- object:

  DGrowthR object.

- data:

  A string indicating the source data to embed. At the moment only
  embeds "od_data"

- ...:

  Additional arguments to pass to the optns argument of uwot

## Value

The result from umap function from uwot package

## See also

[`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html)
