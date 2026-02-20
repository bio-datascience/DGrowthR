# Read individual growth plate

This function receives a file name that contains plate reader OD data.
It returns a data.frame with the parsed data.

## Usage

``` r
read_individual_plate(plate_file, path_target, time_vector = c(), sep = "\t")
```

## Arguments

- plate_file:

  A character variable that is the file name for the plate OD data.
  Should contain a Time column and a column for each well.

- path_target:

  The path to plate_file.

- time_vector:

  A numeric vector that indicates the Time at which each measurement was
  taken.

- sep:

  The delimiter character that separates columns in plate_file.

## Value

A data.frame with the data from plate_file. The Time column contains the
values in time_vector, and the wells are renamed to be plateName_well.
