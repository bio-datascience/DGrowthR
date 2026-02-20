# Instiantiate a DGrowthR object.

This function receives a dataframe with the optical density
measurements, and, optionally and metadata data frame.

## Usage

``` r
DGrowthRFromData(od_data, metadata = NULL, verbose = TRUE)
```

## Arguments

- od_data:

  A data.frame containing all of the optical density data. It should
  have a "Time" column with each row being a timepoint and the rest of
  the columns being an individual growth curve. The names of the growth
  curves should match those in the "curve_id" field in metadata.

- metadata:

  Optional. A data.frame containing a "curve_id" field and all
  additional metadata. If none is provided, then one is created
  automatically containing only the wells from the curve_ids in od_data.

- verbose:

  A logical variable indicating if DGrowthR object should be verbose.

## Value

A DGrowthR object with filled raw_od, od_data, and metadata slots.
