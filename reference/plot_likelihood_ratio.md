# Plot Likelihood Ratios

This function plots a histogram of permuted Likelihood Ratios with the
actual Likelihood Ratio indicated by a vertical line.

## Usage

``` r
plot_likelihood_ratio(object, num_bins = 50)

# S4 method for class 'DGrowthR'
plot_likelihood_ratio(object, num_bins = 50)
```

## Arguments

- object:

  An object of class "DGrowthR".

- num_bins:

  Number of bins to be used in the histogram.

## Value

A ggplot object representing the histogram plot.
