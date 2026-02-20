# Gamma approximated pvalues

A wrapper around the `compute_p_values` function from the permAprox
package using a gamma control distribution.

## Usage

``` r
gamma_pvalues(
  result_list,
  alternative = "greater",
  nullValue = 0,
  fit_thresh = 0.2
)
```

## Arguments

- result_list:

  A list object that is the result of the `multiple_comparisons`
  function with `save_perm_stats == TRUE`.

- alternative:

  A string indicating what the alternative hypothesis should be. Default
  is `greater`.

- nullValue:

  Numeric or character. Specifies the value around which the null
  distribution is centered. If set to "mean" or "median", the per-row
  mean or median of perm_stats is used instead. This allows testing
  against a null hypothesis other than zero or centering based on the
  empirical distribution.

- fit_thresh:

  Numeric. Threshold on empirical p-values below which parametric
  fitting is applied. Default: 0.2 (parametric approximation if the
  empirical p-value is smaller than 0.2).

## Value

Provides a data.frame with the updated the results.
