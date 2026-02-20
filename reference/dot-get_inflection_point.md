# Calculate the inflection point and its tangent of a curve

This function calculates the inflection point (the point of maximum
slope) and the tangent line at this point of a curve. The curve is
defined by the provided dataframe which contains timepoints and
corresponding optical density values.

## Usage

``` r
.get_inflection_point(XX)
```

## Arguments

- XX:

  A dataframe containing the timepoints in a column named "timepoint"
  and the corresponding optical density values in a column named "od".

## Value

A list containing the inflection point and the tangent line data. The
inflection point is returned as a dataframe row from the input
dataframe, with an additional column "curve_id" labeled as
"inflection_point". The tangent line data is a dataframe with timepoints
and corresponding optical density values on the tangent line, with a
column "curve_id" labeled as "inflection_tangent".
