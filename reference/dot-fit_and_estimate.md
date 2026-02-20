# Fit GP and gather growth parameters

Fit GP and gather growth parameters

## Usage

``` r
.fit_and_estimate(
  object,
  model_covariate,
  predict_n_steps,
  downsample_every_n_timepoints,
  sample_posterior_gpfit,
  od_auc_at_t,
  sample_n_curves,
  save_gp_data,
  n_cores,
  detect_polyauxic,
  polyauxic_ratio_threshold,
  polyauxic_parameter_criteria
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
  the second derivative of the mean GP fit. If TRUE, then and additional
  dataframe with the timepoints of the beginning and end of each phase
  is included. If set to TRUE, then save_gp_data must also be set to
  TRUE.

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

Updated DGrowthR object

## See also

[`laGP::predGPsep()`](https://rdrr.io/pkg/laGP/man/predGP.html)
