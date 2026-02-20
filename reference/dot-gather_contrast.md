# Gather the data necesarry to perform a growth comparions.

This function gathers the relevant data for the contrast indicated by
the user. It is an internal function that is used as part of the
growth_comparison method for the DGrowthR class.

## Usage

``` r
.gather_contrast(object, contrast_vector, downsample_every_n_timepoints = 1)
```

## Arguments

- object:

  A DGrowthR object

- contrast_vector:

  A 'vector' object containing three values: 1) The variable on which
  the contrast is made, 2) alternative condition, and 3) reference
  condition

- downsample_every_n_timepoints:

  A numeric value indicating that the OD from every n timepoint should
  be used for GP fit. Might speed up fitting.

## Value

The input object, with the appropiate OD data.

The input object, with the appropiate OD data.
