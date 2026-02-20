# Determine auc

Determine auc

## Usage

``` r
.determine_auc(df, tp_limit = NULL)
```

## Arguments

- df:

  A data.frame with timepoint and mean fields.

- tp_limit:

  A numeric value indicating a timepoint t until which the AUC is
  calculated

## Value

A data.frame with the empirical AUC from the mean GP fit
