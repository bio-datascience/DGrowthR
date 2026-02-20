# Clustering Function

This function performs clustering using the specified algorithm (DBSCAN
or GMM) on UMAP or fPCA coordinates.

## Usage

``` r
clustering(object, ...)

# S4 method for class 'DGrowthR'
clustering(
  object,
  embedding = "umap",
  algorithm = "dbscan",
  k = 0.01,
  eps = NULL
)
```

## Arguments

- object:

  A DGrowthR object containing preprocessed data including UMAP
  coordinates.

- ...:

  Additional parameters to clustering algorithm

- embedding:

  A string indicating if the coordinates from "umap" or "fpca" should be
  used for clustering.

- algorithm:

  A string specifying the clustering algorithm to use ("dbscan" or
  "gmm").

- k:

  Numeric value indicating the number of nearest neighbors (for DBSCAN)
  or number of mixture components (for GMM). If k \< 1, it is
  interpreted as a proportion of the total number of curves.

- eps:

  Numeric value indicating the epsilon distance to use for DBSCAN
  clustering. If NULL, automatic selection based on k is performed.

## Value

Updated DGrowthR object with updated metadata dataframe indicating
cluster membership.

## See also

[`dbscan::dbscan()`](https://rdrr.io/pkg/dbscan/man/dbscan.html),
[`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html)
