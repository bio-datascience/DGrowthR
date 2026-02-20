# Estimate growth parameters

This function estimates a collection of growth parameters from the mean
GP fits. If necessary, it also runs the GP fits.

## Usage

``` r
estimate_growth_parameters(
  object,
  model_covariate = "curve_id",
  predict_n_steps = 100,
  downsample_every_n_timepoints = 1,
  sample_posterior_gpfit = FALSE,
  od_auc_at_t = NULL,
  sample_n_curves = 100,
  save_gp_data = FALSE,
  n_cores = 1,
  detect_polyauxic = FALSE,
  polyauxic_ratio_threshold = 0.2,
  polyauxic_parameter_criteria = "carrying_capacity"
)

# S4 method for class 'DGrowthR'
estimate_growth_parameters(
  object,
  model_covariate = "curve_id",
  predict_n_steps = 100,
  downsample_every_n_timepoints = 1,
  sample_posterior_gpfit = FALSE,
  od_auc_at_t = NULL,
  sample_n_curves = 100,
  save_gp_data = FALSE,
  n_cores = 1,
  detect_polyauxic = FALSE,
  polyauxic_ratio_threshold = 0.2,
  polyauxic_parameter_criteria = "carrying_capacity"
)
```

## Arguments

- object:

  A DGrowthR object containing preprocessed data.

- model_covariate:

  A string indicating a covariate in metadata to pool growth curves for
  GP modelling.

- predict_n_steps:

  A numeric value indicating the number of timepoints to sample from the
  mean posterior of the GP fit.

- downsample_every_n_timepoints:

  A numeric value indicating that the OD from every n timepoint should
  be used for GP fit. Might seep up fitting.

- sample_posterior_gpfit:

  A logical value. Indicates if the posterior GP should be sampled for
  new growth curves and parameters are estimated from these sampled
  curves. Allows to estimate mean and standard deviation of growth
  parameters. Really only makes sense if model_covariate pools more than
  one growth curve per group.

- od_auc_at_t:

  A numeric value indicating a specific timepoint for which the
  predicted OD and AUC should be returned.

- sample_n_curves:

  A numeric values. If sample_posterior_gpfit is TRUE, then
  sample_n_curves are sampled from posterior.

- save_gp_data:

  A logical value, indicating if the mean GP values and GP fit
  parameters should be saved to object.

- n_cores:

  A numeric values indicating the number of cores the user wants to use
  to model curves in parallel.

- detect_polyauxic:

  A logical value indicating if the function should attempt to detect
  polyauxic growth curves. This is done by looking for multiple peaks in
  the second derivative of the mean GP fit. If TRUE, then additional
  columns to the growth parameter table will be added with the relevant
  information.

- polyauxic_ratio_threshold:

  A numeric value. Potential growth phases are assessed by comparing
  their respective carrying capacity. If the carrying capacity of a
  potential second phase is at least this fraction of the carrying
  capacity of the first phase, then the second phase is considered a
  true growth phase. Only relevant if detect_polyauxic is TRUE.

- polyauxic_parameter_criteria:

  A character value. Either "growth_rate" or "carrying_capacity".
  Indicates which parameter should be used to assess the relevance of
  potential growth phases. Only relevant if detect_polyauxic is TRUE.

## Value

Updated DGrowthR object with updated growth_parameters and, if
requested, gpfit_info slots.
