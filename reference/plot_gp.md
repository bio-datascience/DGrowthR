# Plot Gaussian Process

This function plots a specific gaussian process with the 95% confidence
interval.

## Usage

``` r
plot_gp(object, gpfit_id = NULL, show_input_od = FALSE)

# S4 method for class 'DGrowthR'
plot_gp(object, gpfit_id = NULL, show_input_od = FALSE)
```

## Arguments

- object:

  An object of class "DGrowthR"

- gpfit_id:

  The name of the GP model that should be plotted. Should be in the
  gpfit_info\$gpfit_data

- show_input_od:

  A logical value indicating if the od curves used to create the GP
  should be plotted

## Value

A ggplot object showing lineplots of growth curves and GP mean.
