# Predict the values of the fitted model to an object.

This function fits a Gaussian Process model to a DGrowthR object and
then generates predictions from the fitted model. The function supports
"alternative" and "null" models. The "alternative" model includes both
time and contrast as predictors, while the "null" model only includes
time as a predictor.

## Usage

``` r
growth_comparison(
  object,
  comparison_info,
  predict_n_steps = 100,
  downsample_every_n_timepoints = 1,
  save_gp_data = FALSE,
  permutation_test = FALSE,
  n_permutations = NULL,
  n_cores = 1
)

# S4 method for class 'DGrowthR'
growth_comparison(
  object,
  comparison_info,
  predict_n_steps = 100,
  downsample_every_n_timepoints = 1,
  save_gp_data = FALSE,
  permutation_test = FALSE,
  n_permutations = NULL,
  n_cores = 1
)
```

## Arguments

- object:

  A DGrowthR object to which the model should be fitted and then
  predictions should be generated.

- comparison_info:

  character vector describing the contrast to be made. The first value
  indicates the column in metadata where the contrast takes place. The
  third value is taken as the reference condition.

- predict_n_steps:

  A numeric value indicating the number of timepoints to make a
  prediction for with the fitted GP model. More points increases
  estimate accuracy, but adds compute time.

- downsample_every_n_timepoints:

  A numeric value indicating that the OD from every n timepoint should
  be used for GP fit. Might speed up fitting.

- save_gp_data:

  A logical indicating whether the fitted GP models should be stored for
  future analysis.

- permutation_test:

  A logical value indicating if a permutation test should be performed.

- n_permutations:

  A numerical values indicating the number of permutations to build in
  order to gather the null distribution of test statistics.

- n_cores:

  A numeric value indicating the number cores to be used. Useful when
  performing a permutation test.

## Value

An object with the updated growth_comparison slot.
