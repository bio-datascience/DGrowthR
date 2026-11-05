# DGrowthR Development

This checkout is the bio-datascience DGrowthR codebase on the
`codex/dgrowthr-development` branch. The package version associated with the
manuscript under review must remain unchanged. Do not push, merge, or open a
pull request without the user's explicit request; in particular, never push
directly to the bio-datascience main branch.

Read `docs/development_release.md` before release or Git work. Keep GP fitting,
prediction, growth-parameter estimation, and fit QC separate. Preserve the
existing laGP behavior while validating the exact Gaussian backend against
real package curves. A subsampled time grid is selected once for an inference
run and reused for observed and permuted fits. Do not silently resample after
permuting labels.

The Bachelor student's spline work is a separate modelling track; read
`docs/spline_contribution.md` before changing its API. Do not put splines into
`gp_backend` or describe spline uncertainty as a GP credible interval.

Before calling development work complete, run focused tests and a package
check when dependencies are available. Report missing dependencies and failed
checks explicitly. Preserve unrelated changes in this shared checkout.
