# Preprocess Growth Data

This function reads growth data from a file, applies several optional
preprocessing steps, and saves the processed data.

## Usage

``` r
preprocess_data(
  object,
  baseline = 1,
  log_transform = TRUE,
  skip_first_n_timepoints = 0
)

# S4 method for class 'DGrowthR'
preprocess_data(
  object,
  baseline = 1,
  log_transform = TRUE,
  skip_first_n_timepoints = 0
)
```

## Arguments

- object:

  A DGrowthR object containing the growth data.

- baseline:

  A numeric value specifying timepoint that will be used for baseline
  correction. If a vector is provided, then the minimum OD of the
  provided timepoints is used as baseline.

- log_transform:

  A logical value indicating whether to perform log transformation on
  optical density measurements. This enables accurate estimation of
  growth parameters.

- skip_first_n_timepoints:

  A numeric value the number of timepoints to skip for all growth
  curves. The analysis would start at skip_first_n_timepoints + 1.

## Value

The processed DGrowthR object with updated od_data.
