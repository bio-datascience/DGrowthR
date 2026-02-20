# Detect polyauxic growth curves This function is heavily based on the implementation available on AMiGA.

Detect polyauxic growth curves This function is heavily based on the
implementation available on AMiGA.

## Usage

``` r
.detect_polyauxic(df, thresh, varb_detect)
```

## Arguments

- df:

  A data.frame with timepoint, mean, first and second derivative fields.

- thresh:

  A numeric value indicating the minimum fraction of the carrying
  capacity of the first phase that the second phase must have to be
  considered a true growth phase.

- varb_detect:

  A character value. Either "growth_rate" or "carrying_capacity".
  Indicates which parameter should be used to assess the relevance of
  potential growth phases. Only relevant if detect_polyauxic is TRUE.

## Value

A data frame.
