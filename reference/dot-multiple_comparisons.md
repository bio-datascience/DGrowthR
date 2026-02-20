# A wrapper for performing multiple comparisons that can catch errors.

A wrapper for performing multiple comparisons that can catch errors.

## Usage

``` r
.multiple_comparisons(
  index,
  comparison_list,
  object,
  predict_n_steps = 100,
  downsample_every_n_timepoints = 1,
  save_gp_data = FALSE,
  permutation_test = FALSE,
  n_permutations = NULL,
  n_cores = 1,
  save_perm_stats = FALSE,
  write_to_tsv = FALSE,
  filename = NULL
)
```

## Arguments

- comparison_list:

  A list of contrasts.

- object:

  A DGrowthR object to which the model should be fitted and then
  predictions should be generated.

- predict_n_steps:

  A logical

- downsample_every_n_timepoints:

  A numeric value indicating that the OD from every n timepoint should
  be used for GP fit. Might seep up fitting.

- save_gp_data:

  A logical value, indicating if the mean GP values and GP fit
  parameters should be saved to object.

- permutation_test:

  A logical value indicating if a permutation test should be performed.

- n_permutations:

  A numerical values indicating the number of permutations to build in
  order to gather the null distribution of test statistics.
