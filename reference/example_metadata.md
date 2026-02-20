# Example DGrowthR Metadata

A dataset containing metadata information related to different
treatments and their properties.

## Usage

``` r
example_metadata
```

## Format

A data frame with 150 rows and 5 variables:

- curve_id:

  A character string indicating the ID or name of the curve.

- Treatment:

  A character string indicating the type and, the concentration of the
  treatment.

- Conc.:

  A numeric value indicating the concentration of the treatment.

- Product:

  A character string describing the type of product or compound related
  to the treatment. For example, "Cancer Medication", "Water", "DMSO",
  or "HIV Medication".

- abx:

  A character string indicating whether an antibiotic (or other relevant
  substance) is present. Values are either "abx" or "no_abx".

## Examples

``` r
head(DGrowthR::example_metadata)
#>          curve_id   Treatment Conc.           Product abx
#> 1   Cancer_0.05_1 Cancer_0.05  0.05 Cancer Medication abx
#> 1.1 Cancer_0.05_2 Cancer_0.05  0.05 Cancer Medication abx
#> 1.2 Cancer_0.05_3 Cancer_0.05  0.05 Cancer Medication abx
#> 1.3 Cancer_0.05_4 Cancer_0.05  0.05 Cancer Medication abx
#> 1.4 Cancer_0.05_5 Cancer_0.05  0.05 Cancer Medication abx
#> 1.5 Cancer_0.05_6 Cancer_0.05  0.05 Cancer Medication abx
```
