# Spline modelling contribution

This is a scoped starting point for the Bachelor student contributing to the
DGrowthR development branch. The existing laGP and exact Gaussian GP behavior
is the reference baseline; spline work must not change its numerical results.

## Statistical contract to decide first

- State whether the spline estimates each replicate separately or a pooled
  group-mean curve. Replicates are biological measurements, not extra time
  coordinates; do not silently average or weight them without documenting the
  estimator.
- State the basis, knot policy, smoothing/regularization parameter selection,
  boundary behavior, and handling of irregular timepoints and flat curves.
- Define first and second derivatives analytically from the fitted spline.
- Distinguish a fitted mean curve from uncertainty about the mean and from a
  predictive interval for observations. Do not label a residual band as a GP
  credible interval.
- Decide whether the fitting objective yields a likelihood suitable for the
  existing likelihood-ratio permutation test. If not, use a separate test
  statistic with its own null permutation distribution and documentation.

## Code boundary

The current `fit_gp()`, `predict_gp()`, `fit_predict_gp()`, and `assess_gp_fit()`
are GP-specific. Do not add `"spline"` to `gp_backend` or feed a spline fit into
the GP prediction or QC functions. A spline implementation can first expose
internal fit/predict functions that return a serializable fit and a prediction
table with `timepoint`, `mean`, `first_derivative`, and `second_derivative`.
Once validated, add a model-neutral dispatch layer or explicit spline API;
keep the GP API as a stable compatibility wrapper.

## Initial validation

1. Fit real package curves with 15-150 timepoints, including 3-6 replicates,
   flat curves, and irregular time grids.
2. Compare fitted means and derivatives with exact GP and finite differences,
   especially near the start/end and at multiple steep ascents.
3. Verify growth-parameter estimates, including the maximum growth rate, on a
   dense common prediction grid.
4. Report runtime and memory separately for fitting and prediction.
5. Exercise subsampled inference on one fixed selected time grid; never choose
   a different grid after permuting labels.

The first deliverable should be a small design note and a tested prototype,
not a change to the public GP backend switch.
