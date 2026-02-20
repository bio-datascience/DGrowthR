# Estimate lag phase

Estimate lag phase

## Usage

``` r
.determine_lag_phase(df, max_gr_data)
```

## Arguments

- df:

  A data.frame with timepoint and second_derivative fields.

- max_gr_data:

  A data.frame with the output of .determine_max_growth_rate().

## Value

A data.frame with the timepoint at the end of lag-phase
