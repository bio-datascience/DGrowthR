# Plot Optical Density

This function plots the raw optical density curves, with the possibility
of specifying color and facet variables

## Usage

``` r
plot_growth_curves(object, color = NULL, facet = NULL, data = "od_data")

# S4 method for class 'DGrowthR'
plot_growth_curves(object, color = NULL, facet = NULL, data = "od_data")
```

## Arguments

- object:

  An object of class "DGrowthR"

- color:

  A column of metadata to color growth curves

- facet:

  A column of metadata to facet growth curves plots

- data:

  A character value indicating which data should be plotted. If "raw_od"
  then unprocessed data is plotted. If "od_data" (default) then the
  pre-processed data is plotted if present.

## Value

A ggplot object showing lineplots of growth curves.
