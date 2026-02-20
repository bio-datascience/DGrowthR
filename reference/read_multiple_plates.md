# Read a collection growth plate data.

This function receives a path to a folder that contains data for several
plates. Each plate should have its own file.

## Usage

``` r
read_multiple_plates(path_target, time_vector = c(), sep = "\t")
```

## Arguments

- path_target:

  The path to plate_file.

- time_vector:

  A numeric vector that indicates the Time at which each measurement was
  taken. If empty then the Time column from one of the plates is taken
  as reference.

- sep:

  The delimiter character that separates columns of each file in
  path_target

## Value

A list of two dataframes. An "od_data" dataframe with a common time
column and each growth curve in a separate column named with it's
"curve_id" (plateName_well). Also a "metadata" dataframe indicating the
well and the original plate from which each "curve_id" is gathered from.
