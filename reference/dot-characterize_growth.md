# Estimate growth parameters

Estimate growth parameters

## Usage

``` r
.characterize_growth(
  df,
  od_auc_at_t = NULL,
  detect_polyauxic = FALSE,
  polyauxic_ratio_threshold = 0.2,
  polyauxic_parameter_criteria = "carrying_capacity"
)
```

## Arguments

- df:

  A data.frame with timepoint and mean fields.

- od_auc_at_t:

  A numeric value for which od and AUC are calculated for

- detect_polyauxic:

  A logical value indicating if the function should attempt to detect
  polyauxic growth curves. This is done by looking for multiple peaks in
  the second derivative of the mean GP fit. If TRUE, then and additional
  dataframe with the timepoints of the beginning and end of each phase
  is included.

- polyauxic_ratio_threshold:

  A numeric value. Potential growth phases are assessed by comparing
  their respective carrying capacity. If the carrying capacity of a
  potential second phase is at least this fraction

- polyauxic_parameter_criteria:

  A character value. Either "growth_rate" or "carrying_capacity".
  Indicates which parameter should be used to assess the relevance of
  potential growth phases. Only relevant if detect_polyauxic is TRUE.

## Value

A data.frame with the requested growth parameters for the curve given as
input
