# Gather specific growth data from the DGrowthR object

Extracts the requested growth curves according to the requested metadata
fields and variables.

## Usage

``` r
gather_od_data(
  object,
  metadata_field,
  field_value,
  downsample_every_n_timepoints = 1
)

# S4 method for class 'DGrowthR'
gather_od_data(
  object,
  metadata_field,
  field_value,
  downsample_every_n_timepoints = 1
)
```

## Arguments

- object:

  An object of class "DGrowthR".

- metadata_field:

  A character string specifying one of the fields in the metadata from
  which the variables will be searched in.

- field_value:

  The value of the metadata_field that for which the relevant growth
  curves will be extracted.

- downsample_every_n_timepoints:

  A numeric value indicating that the OD from every n timepoint should
  be used for GP fit. Might seep up fitting.

## Value

A data.frame containing the optical density data requested.
