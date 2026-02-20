# Determine optimal value of epsilon with elbow-detection

This function determines the optimal value of epsilon for density-based
clustering based on the sorted K-nearest neighbor distance curve. The
optimal epsilon is the 'elbow' of the curve and is found as the point
with the largest distance to the (imaginary) line defined by the first
and last point.

## Usage

``` r
.determine_optimal_eps(X, k)
```

## Arguments

- X:

  numeric matrix of coordinates in UMAP space

- k:

  numeric value that indicates the K nearest neighoburs to be used

## Value

numeric value corresponding to the optimal eps value
