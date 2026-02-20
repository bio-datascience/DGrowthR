# Perform a permutation test for the DGrowthR object

This function performs a permutation test for the DGrowthR object. It
shuffles the "contrast" column in the OD data and calculates the
difference in log-likelihoods between the alternative and null models
for each permutation. The permutation test helps assess the significance
of the difference in models.

## Usage

``` r
perm_test(
  object,
  num_perms = 1000,
  n.cores = 1,
  predict_n_steps = 100,
  gp.delete = TRUE,
  ...
)

# S4 method for class 'DGrowthR'
perm_test(
  object,
  num_perms = 1000,
  n.cores = 1,
  predict_n_steps = 100,
  gp.delete = TRUE,
  ...
)
```

## Arguments

- object:

  A DGrowthR object which contains the Gaussian Process regression
  results.

- num_perms:

  The number of permutations to perform. Defaults to 1000.

- n.cores:

  The number of CPU cores to use for parallel processing. Defaults to 1
  (no parallel processing).

- predict_n_steps:

  A numeric value indicating the number of timepoints to make a
  prediction for with the fitted GP model. More points increases
  estimate accuracy, but adds compute time.

- gp.delete:

  A logical value indicating whether to delete the GPsep identifiers
  after the permutation test. Defaults to TRUE.

- ...:

  Additional arguments to be passed to the GP modelling function.

## Value

A modified DGrowthR object with the calculated log-likelihood
differences for each permutation.
