# Future DGrowthR Speed Optimizations

This note records speed strategies that should preserve DGrowthR's numerical
results. The intent is to survive thread summaries and keep future work focused
on avoiding repeated work rather than changing the statistical model.

## Current State

- DGrowthR now has an exact Gaussian GP backend with analytic first and second
  posterior-mean derivatives.
- `subsample_timepoints()` defaults to a Krause/Singh/Guestrin-style greedy GP
  mutual-information criterion for selecting inference timepoints.
- The GP itself does not require uniformly spaced timepoints.
- Typical trajectories are short, often 15 timepoints and usually no more than
  about 150 timepoints.
- Replicate curves are usually measured on the same time grid.

## Optimization Principles

1. Preserve the GP and numerical results whenever possible.
2. Prefer caching and reuse over approximate math.
3. Treat repeated time grids as a first-class object.
4. Do not compromise permutation-test fairness: all permutations and comparison
   groups should see the same inference grid when subsampling is used.
5. Avoid repeated fitting where the design matrix, response, and hyperparameter
   policy are identical.

## Priority Ideas

### 0. Possible Fourier GP Backend From EFGP

Greengard, Rachh, and Barnett's EFGP method is a different way to compute GP
posterior means, not a subsampling method. The reference MATLAB code is:

```text
https://github.com/flatironinstitute/gp-shootout
```

The DGrowthR-relevant files are:

```text
algs/EFGP/EFGP.m
algs/EFGP/efgp1d.m
algs/EFGP/get_xis.m
kernels/SE_ker.m
```

The method approximates a stationary kernel by equispaced Fourier features:

```text
k(t - t') ~= sum_j h * khat(hj) * exp(2*pi*i*h*j*(t - t'))
```

Then it solves the weight-space system:

```text
(Phi* Phi + noise I) beta = Phi* y
```

instead of the dense function-space system:

```text
(K + noise I) alpha = y
```

For DGrowthR's exact Gaussian parameterization:

```text
k(t, t') = exp(-(t - t')^2 / d)
```

this corresponds to the paper's squared-exponential length scale:

```text
ell = sqrt(d / 2)
```

The normalized exact GP mean can use:

```text
prior variance = 1
noise variance = relative nugget
```

because the profiled process variance cancels out of the posterior mean when
both prior covariance and noise are scaled by the same `sigma2`.

Important caveats:

- The paper/code target fast posterior means for fixed hyperparameters.
- Hyperparameter likelihoods, determinants, and likelihood gradients are not
  implemented in the EFGP paper; the paper explicitly leaves this for future
  work.
- DGrowthR still needs exact dense GP fitting for lengthscale/nugget
  optimization unless we develop an approximate Fourier likelihood backend.
- For typical DGrowthR curves with 15-150 timepoints, dense exact Cholesky may
  still be faster than EFGP because EFGP has FFT/NUFFT/CG overhead.
- EFGP could be useful for large pooled raw-replicate fits, full screens,
  repeated prediction on many target grids, or as a fixed-hyperparameter
  approximate backend.

Practical implementation route:

1. Add `gp_backend = "fourier_gaussian"` only for fixed hyperparameters first.
2. Implement a pure-R dense Fourier-feature prototype for 1D, without NUFFT,
   because DGrowthR time grids are small.
3. Compare posterior mean and analytic derivatives to `exact_gaussian`.
4. Only consider C/FFT/FINUFFT-style acceleration if the prototype wins for
   realistic DGrowthR workloads.

### 1. Cache Time-Grid Geometry

If all curves share the same timepoints, then squared-distance matrices are
identical across fits:

```text
D_ij = (t_i - t_j)^2
```

Cache and reuse:

- fitted-grid squared distances,
- prediction-grid squared distances,
- cross-distances between fitted and prediction grids,
- selected inference grids from `subsample_timepoints()`.

This preserves results exactly and removes repeated `O(n^2)` setup costs.

### 2. Reuse the Inference Subsampling Grid

The default information-selection strategy depends only on the time grid and
kernel hyperparameters, not OD values. Compute the selected timepoints once for
an inference run, then reuse the same thinned DGrowthR object for all
`growth_comparison()` and permutation calls.

This is important for unbiased/fair permutation strategies because every label
permutation should use the same observed time grid.

### 3. Aggregate Replicates When Estimating Group Mean Curves

When replicates share identical timepoints, averaging replicates per
group/timepoint is an unbiased estimate of the group mean trajectory:

```text
\bar y_g(t) = (1 / r_g) sum_{j=1}^{r_g} y_{gj}(t)
```

This reduces the GP matrix size from roughly:

```text
n_timepoints * n_replicates
```

to:

```text
n_timepoints
```

Because exact GP fitting is cubic in matrix size, this can be a major speed win.
Use this only when the estimand is the group mean trajectory, not individual
replicate trajectories.

### 4. Memoize Repeated GP Fits

Create a fit cache keyed by:

```text
backend
time grid
prediction grid
response vector hash
contrast design
hyperparameter policy
lengthscale/nugget settings
replicate strategy
```

This should be most useful inside `perm_test()`, where many null/alternative
fits can share identical inputs or geometry.

### 5. Avoid Full Posterior Covariance Unless Needed

Use mean-only predictions by default. Full covariance is only needed for
posterior sampling.

Keep:

```r
complete_sigma = FALSE
```

unless `sample_posterior_gpfit = TRUE` or another workflow explicitly needs
posterior draws.

### 6. Use Analytic Exact-GP Derivatives Everywhere Possible

The exact Gaussian backend has analytic first and second derivatives. Prefer
these over finite-difference derivatives and avoid increasing prediction-grid
resolution solely to stabilize finite differences.

### 7. Reduce R Loop Overhead in Permutation Workflows

Keep the same objective and optimizer, but move repeated permutation bookkeeping
toward vectorized matrix code or compiled code. This should come after caching,
because repeated unnecessary GP fits are likely the larger cost.

## Suggested Implementation Order

1. Make selected inference timepoints easy to extract and reuse from
   `subsample_timepoints()`.
2. Add an inference helper that creates one thinned DGrowthR object and records
   the selected grid.
3. Add explicit time-grid distance caches for exact GP fitting/prediction.
4. Add memoization around exact GP fits inside `growth_comparison()` and
   `perm_test()`.
5. Profile permutation workflows and only then consider compiled permutation
   bookkeeping.

## Reference

Krause, A., Singh, A., & Guestrin, C. (2008). Near-optimal sensor placements in
Gaussian processes: Theory, efficient algorithms and empirical studies. Journal
of Machine Learning Research, 9, 235-284.
