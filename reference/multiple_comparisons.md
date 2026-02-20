# Multiple Comparisons

This function processes multiple contrasts, performs an analysis, and
returns a data frame with the results.

Key features include:

- **Multiple contrast processing**: Ability to process and analyze
  multiple contrasts.

- **Comprehensive Analysis**: Performs a range of analyses based on user
  preferences, including Likelihood Ratio, empirical p-value, area under
  the curve (AUC) for both condition and control, maximum growth rate of
  the condition, maximum difference, and maximum optical density (OD)
  for both condition and control.

- **Output and Display Formats**: Provides an option to save results to
  a TSV file.

## Usage

``` r
multiple_comparisons(
  object,
  comparison_list,
  predict_n_steps = 50,
  downsample_every_n_timepoints = 1,
  permutation_test = FALSE,
  n_permutations = NULL,
  n_cores = 1,
  write_to_tsv = FALSE,
  save_perm_stats = FALSE,
  filename = NULL
)
```

## Arguments

- object:

  DGrowthR object that has been instantiated.

- comparison_list:

  A list of of vectors. Each vector contains the information for an
  individual contrast.

- predict_n_steps:

  A numeric value specifying the number of time points for prediction.
  Defaults to 50.

- downsample_every_n_timepoints:

  A numeric value that indicates how timepoints are sampled to increase
  fit speed

- permutation_test:

  A logical value indicating if a permutation test should be performed.

- n_permutations:

  Number of permutations for the permutation test (default: 100).

- n_cores:

  Number of cores to use for the permutation test (default: 1).

- write_to_tsv:

  A logical value inndicating if the results should be save to a tsv
  file.

- save_perm_stats:

  A logical value indicating whether the permuted test statistics should
  be stored in the DGrowthR object.

- filename:

  Name of the TSV file.

## Value

Provides a data.frame summarizing the results.
