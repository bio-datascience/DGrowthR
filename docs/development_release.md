# DGrowthR development release

The `codex/dgrowthr-development` branch is the integrated development line for
DGrowthR 1.1 (`DESCRIPTION` version `1.1.0.9000`). It is separate from the
package version under manuscript review. Do not push this branch to the main
repository or open upstream pull requests without an explicit decision.

## Scope of the development line

- Exact Gaussian GP with laGP-style initialization, optimized lengthscale and
  nugget, serializable fitted state, symbolic posterior mean, and analytic
  first and second derivatives.
- Separate `fit_gp()`, `predict_gp()`, `estimate_growth_parameters()`, and
  `assess_gp_fit()` steps, while retaining `fit_predict_gp()` for existing code.
- Explicit timepoint subsampling with a Krause-style information criterion as
  the default, plus even and GP-landmark alternatives. One selected inference
  grid is used for observed fits and every permutation. This changes the
  inference data, not the GP model family.
- Replicate-aware raw and aggregate fitting strategies, with exact-GP
  aggregate fitting using heteroskedastic nugget weights.
- A separate student-led spline modelling track. It should share growth-curve
  input and prediction conventions, but remain distinct from Gaussian-process
  fitting and GP-specific uncertainty diagnostics until its statistical
  contract has been defined and validated. See `docs/spline_contribution.md`.

## Release gates before a fork push

1. Install missing R dependencies and run `R CMD check` on a built tarball.
2. Test fit/predict/QC, subsampling, comparison, and permutation paths together
   on real package growth curves, including grouped replicates and irregular
   timepoints.
3. Compare exact-GP and laGP numerical behavior and runtimes on the package's
   representative short trajectories.
4. Review API documentation and compatibility with existing scripts.
5. Confirm the fork remote and push only to the user's fork.

Upstream pull requests can later be split by behavior: exact GP backend,
workflow API and QC, subsampling and inference integration, then a spline
modelling extension after its own validation. This document is a release plan,
not authorization to push or open a PR.
