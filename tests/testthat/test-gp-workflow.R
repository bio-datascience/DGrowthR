test_that("stored fits predict and survive serialization", {
  x <- seq(0, 2, length.out = 12)
  object <- methods::new(
    "DGrowthR",
    od_data = data.frame(curve_id = rep(c("a", "b"), each = 24),
                         timepoint = rep(rep(x, 2), 2),
                         timepoint_n = rep(rep(seq_along(x), 2), 2),
                         od = rep(c(sin(x), sin(x) + 0.01), 2)),
    metadata = data.frame(curve_id = c("a", "b")),
    preprocessed = TRUE, log_od = TRUE)

  for (backend in c("lagp", "exact_gaussian")) {
    fitted <- fit_gp(object, gp_backend = backend,
                     optimize_hyperparameters = FALSE)
    expect_equal(length(fitted@gpfit_info$fits), 2)
    expect_true(all(vapply(fitted@gpfit_info$fits, function(x) is.null(x$fit$gpsepi),
                           logical(1))))
    restored <- unserialize(serialize(fitted, NULL))
    predicted <- predict_gp(restored, prediction_timepoints = x)
    expect_equal(nrow(predicted@gpfit_info$gpfit_data), 24)
    expect_true(all(c("first_derivative", "second_derivative") %in%
                    names(predicted@gpfit_info$gpfit_data)))
    qc <- assess_gp_fit(restored)
    expect_equal(nrow(qc$summary), 2)
    expect_true(all(qc$summary$rmse < 0.1))
    expect_equal(unique(qc$points$interval_type), "predictive")
    expect_equal(nrow(assess_gp_fit(restored, interval_type = "latent")$summary), 2)
    growth <- estimate_growth_parameters(predicted, predict_n_steps = 12)
    expect_equal(nrow(growth@growth_parameters), 2)
    expect_equal(vapply(growth@gpfit_info$fits, function(x) x$fit$lengthscale[1],
                        numeric(1)),
                 vapply(fitted@gpfit_info$fits, function(x) x$fit$lengthscale[1],
                        numeric(1)))
  }
})

test_that("flat groups retain a status without stopping other fits", {
  x <- seq(0, 1, length.out = 8)
  object <- methods::new(
    "DGrowthR",
    od_data = data.frame(curve_id = rep(c("flat", "growth"), each = 8),
                         timepoint = rep(x, 2),
                         timepoint_n = rep(seq_along(x), 2),
                         od = c(rep(1, 8), sin(x))),
    metadata = data.frame(curve_id = c("flat", "growth")))
  fitted <- fit_gp(object, gp_backend = "exact_gaussian",
                   optimize_hyperparameters = FALSE)
  expect_equal(fitted@gpfit_info$fits$flat$status, "flat_curve")
  expect_equal(fitted@gpfit_info$fits$growth$status, "fitted")
  expect_equal(nrow(assess_gp_fit(fitted)$summary), 2)
})

test_that("QC distinguishes isolated and sustained deviations on irregular times", {
  x <- c(0, 0.1, 0.2, 0.7, 1.3, 2.0, 2.8, 3.7)
  object <- methods::new(
    "DGrowthR",
    od_data = data.frame(curve_id = "a", timepoint = x,
                         timepoint_n = seq_along(x), od = sin(x)),
    metadata = data.frame(curve_id = "a"))
  fitted <- fit_gp(object, gp_backend = "exact_gaussian",
                   optimize_hyperparameters = FALSE)
  isolated <- fitted
  isolated@gpfit_info$fits$a$observations$od[4] <- 10
  isolated_qc <- assess_gp_fit(isolated)
  expect_equal(isolated_qc$summary$status, "pass")
  expect_equal(isolated_qc$runs$n_points, 1)

  sustained <- fitted
  sustained@gpfit_info$fits$a$observations$od[4:6] <- 10
  sustained_qc <- assess_gp_fit(sustained)
  expect_equal(sustained_qc$summary$status, "fail")
  expect_equal(sustained_qc$summary$longest_outside_duration, 1.3)
  expect_equal(sustained_qc$summary$longest_outside_direction, "above")
})

test_that("subsampling is shared with fitting and cannot be applied twice", {
  x <- seq(0, 2, length.out = 12)
  object <- methods::new(
    "DGrowthR",
    od_data = data.frame(curve_id = "a", timepoint = x,
                         timepoint_n = seq_along(x), od = sin(x)),
    metadata = data.frame(curve_id = "a"))
  thinned <- subsample_timepoints(object, n_timepoints = 6, method = "information")
  selected <- attr(thinned@od_data, "subsampled_timepoints")
  fitted <- fit_gp(thinned, gp_backend = "exact_gaussian",
                   optimize_hyperparameters = FALSE)
  expect_equal(fitted@gpfit_info$gpfit_parameters$subsampled_timepoints, selected)
  expect_equal(sort(unique(fitted@gpfit_info$fits$a$observations$timepoint)), selected)
  expect_error(fit_gp(thinned, downsample_every_n_timepoints = 2),
               "already subsampled")
})
