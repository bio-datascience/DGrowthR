# Pad inflection points and directions This function is heavily based on the implementation available on AMiGA.

Pad inflection points and directions This function is heavily based on
the implementation available on AMiGA.

## Usage

``` r
.pad_inflections(ips, its, edge = 1, len = NULL)
```

## Arguments

- ips:

  A numeric vector of inflection points (indices).

- its:

  A numeric vector of inflection directions (1 for positive, -1 for
  negative).

- edge:

  A numeric value indicating which edge to pad (1 for left, -1 for
  right).

- len:

  A numeric value indicating the total length of the data vector (used
  for right edge padding).

## Value

A list containing the padded inflection points and directions.
