# This S4 class represents a DGrowthR model.

It is used to store and validate data and parameters necessary for
DGrowthR analysis.

## Slots

- `od_data`:

  A data.frame expected to contain the following columns: "timepoint",
  "od", and "curve_id". These represent the time points, and optical
  density (od), and an identifier for each curve respectively.

- `metadata`:

  A data.frame that contains the "curve_id" column along with any
  metadata variables

- `raw_od`:

  Optical density information before pre-processing or re-formatting.

- `verbose`:

  A logical value indicating whether the output should be suppressed.
  Defaults to TRUE, i.e., output is not suppressed.

- `fpca`:

  An fpca object as the output of FPCA from fdapace

- `umap_coord`:

  Coordinates for UMAP embedding

- `cluster_assignment`:

  A data.frame matching curve_ids to the cluster assigned by
  clustering() function.

- `gpfit_info`:

  A list with the estimated mean gp process and the parameters used for
  fitting.

- `growth_comparison`:

  The results of a growth comparison.

- `log_od`:

  A logical indicating if the OD data has been logged.

- `preprocessed`:

  A logical variable indicating if the od data has been preprocessed
  with the preprocess_data() function.

- `growth_parameters`:

  A data.frame with the estimated growth parameters after gp modelling

## Validity

For an object to be a valid DGrowthR, it must meet the following
criteria:

- The "od_data" slot must be a data frame containing the columns
  "timepoint", "timepoint_n", "curve_id", and "od".

- The "metadata" slot must be a data frame with the "curve_id" column.

- All of the "curve_id"s in the "od_data" data frame must be present in
  the "curve_id" column of the metadata.
