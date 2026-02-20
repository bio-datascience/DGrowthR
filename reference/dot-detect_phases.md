# Iteratively detect phases This function is heavily based on the implementation available on AMiGA.

Iteratively detect phases This function is heavily based on the
implementation available on AMiGA.

## Usage

``` r
.detect_phases(df, ratio_threshold, varb = "K", second_varb = "r")
```

## Arguments

- df:

  A data.frame with the parameters of each potential growth phase.

- ratio_threshold:

  A numeric value indicating the minimum fraction of the carrying
  capacity of the first phase that the second phase must have to be
  considered a true growth phase.

## Value

A data frame with the detected growth phases.
