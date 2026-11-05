# Either use the provided example od_data
od_data <- DGrowthR::example_od_data
metadata <- DGrowthR::example_metadata
contrast <- c("Treatment", "Cancer_0.05", "DMSO_0.05")

# Create an instance of the DGrowthR class
test_object <- DGrowthRFromData(od_data = od_data, metadata = metadata)

test_that("DGrowthR S4: class can be created with valid inputs", {
 
  DGrowthR_obj <- DGrowthRFromData(od_data = od_data, metadata = metadata)

  expect_s4_class(DGrowthR_obj, "DGrowthR")
  expect_true(is.list(slot(DGrowthR_obj, "growth_comparison")))
  expect_true(is.list(slot(DGrowthR_obj, "gpfit_info")))
  expect_true(is.list(slot(DGrowthR_obj, "growth_comparison")))
  expect_true(is.list(slot(DGrowthR_obj, "fpca")))
  expect_true(is.data.frame(slot(DGrowthR_obj, "umap_coord")))
  expect_false(slot(DGrowthR_obj, "log_od"))
})

test_that("DGrowthR S4: class cannot be created with invalid inputs", {
  # Sample invalid data (missing 'Time' column)
  od_data_invalid <- data.frame(
    sample_1 = c(1, 2, 3),
    sample_2 = c(1, 2, 3)
  )

  expect_error(
    DGrowthRFromData(od_data_invalid, metadata),
    'No column called "Time" in od_data'
  )
})


# Test that we can do a growth comparison
test_that("growth_comparison works", {
  result <- growth_comparison(object=test_object, predict_n_steps = 10, comparison_info = contrast)
  expect_true(is(result, "DGrowthR"))
})

test_that("estimate_growth_parameters works", {
  test_object_pre <- preprocess_data(test_object)
  test_object_params <- estimate_growth_parameters(test_object_pre, save_gp_data=FALSE)
  expect_true(nrow(test_object_params@growth_parameters) == nrow(metadata))
})

test_that("estimate_growth_parameters by group works", {
  test_object_pre <- preprocess_data(test_object)
  test_object_params <- estimate_growth_parameters(test_object_pre, model_covariate="Treatment", predict_n_steps = 10)
  expect_true(nrow(test_object_params@growth_parameters) == 5)
})


test_that("estimate_growth_parameters with posterior sampling works", {
  test_object_pre <- preprocess_data(test_object)
  test_object_params <- estimate_growth_parameters(test_object_pre, 
                                                   model_covariate="Treatment", 
                                                   sample_posterior_gpfit=TRUE, 
                                                   predict_n_steps = 10,
                                                   sample_n_curves=5,
                                                   save_gp_data = FALSE)
  expect_true("max_growth_time.mean" %in% colnames(test_object_params@growth_parameters))
})

test_that("maximum growth rate uses earliest timepoint for exact ties", {
  df <- data.frame(
    timepoint = 1:6,
    mean = c(0, 1, 2, 2.5, 3.5, 4.5),
    first_derivative = c(0.1, 1, 0.2, 1, 0.3, NA),
    second_derivative = c(0.1, 0.2, -0.1, 0.3, -0.2, NA)
  )

  max_growth_rate <- .determine_max_growth_rate(df, .determine_max_growth(df))

  expect_equal(nrow(max_growth_rate), 1)
  expect_equal(max_growth_rate$max_growth_rate, 1)
  expect_equal(max_growth_rate$max_growth_rate_timepoint, 2)
})

test_that("growth characterization returns one row under tied extrema", {
  df <- data.frame(
    timepoint = 1:6,
    mean = c(0, 1, 2, 2, 1, 1),
    first_derivative = c(0.1, 0.3, 0.3, -0.2, -0.2, NA),
    second_derivative = c(0.5, 0.5, -0.1, -0.4, -0.4, NA)
  )

  growth_parameters <- .characterize_growth(df)

  expect_equal(nrow(growth_parameters), 1)
  expect_equal(growth_parameters$max_growth_time, 3)
  expect_equal(growth_parameters$max_death_rate_timepoint, 4)
  expect_equal(growth_parameters$lag_phase_timepoint, 1)
  expect_equal(growth_parameters$stationary_phase_timepoint, 3)
})

test_that("exact Gaussian raw grouped replicates return one row per group", {
  test_object_pre <- preprocess_data(test_object)
  selected_treatments <- c("Cancer_0.05", "DMSO_0.05")
  keep_metadata <- test_object_pre@metadata %>%
    filter(Treatment %in% selected_treatments)

  slot(test_object_pre, "metadata") <- keep_metadata
  slot(test_object_pre, "od_data") <- test_object_pre@od_data %>%
    filter(curve_id %in% keep_metadata$curve_id)

  test_object_params <- estimate_growth_parameters(
    test_object_pre,
    model_covariate = "Treatment",
    predict_n_steps = 10,
    save_gp_data = FALSE,
    n_cores = 1,
    gp_backend = "exact_gaussian",
    gp_replicate_strategy = "raw",
    gp_nugget = 1e-6
  )

  expect_equal(nrow(test_object_params@growth_parameters), length(selected_treatments))
  expect_equal(sort(test_object_params@growth_parameters$gpfit_id), sort(selected_treatments))
})

test_that("estimate_growth_parameters handles a single exact Gaussian raw group", {
  test_object_pre <- preprocess_data(test_object)
  selected_treatment <- "Cancer_0.05"
  keep_metadata <- test_object_pre@metadata %>%
    filter(Treatment == selected_treatment)

  slot(test_object_pre, "metadata") <- keep_metadata
  slot(test_object_pre, "od_data") <- test_object_pre@od_data %>%
    filter(curve_id %in% keep_metadata$curve_id)

  test_object_params <- estimate_growth_parameters(
    test_object_pre,
    model_covariate = "Treatment",
    predict_n_steps = 10,
    save_gp_data = TRUE,
    n_cores = 1,
    gp_backend = "exact_gaussian",
    gp_replicate_strategy = "raw",
    gp_nugget = 1e-6,
    optimize_hyperparameters = FALSE
  )

  expect_s4_class(test_object_params, "DGrowthR")
  expect_true(is.data.frame(test_object_params@growth_parameters))
  expect_equal(nrow(test_object_params@growth_parameters), 1)
  expect_equal(test_object_params@growth_parameters$gpfit_id, selected_treatment)
  expect_true(is.data.frame(test_object_params@gpfit_info$gpfit_data))
  expect_true(length(test_object_params@gpfit_info$gp_symbolic_mean) == 1)
})

test_that("exact Gaussian GP backend returns analytic derivatives and symbolic mean", {
  test_object_pre <- preprocess_data(test_object)
  test_object_params <- estimate_growth_parameters(
    test_object_pre,
    model_covariate = "Treatment",
    predict_n_steps = 10,
    save_gp_data = TRUE,
    gp_backend = "exact_gaussian",
    gp_nugget = 1e-6
  )

  expect_true(nrow(test_object_params@growth_parameters) == 5)
  expect_true("first_derivative" %in% colnames(test_object_params@gpfit_info$gpfit_data))
  expect_true("second_derivative" %in% colnames(test_object_params@gpfit_info$gpfit_data))
  expect_equal(test_object_params@gpfit_info$gpfit_parameters$gp_backend, "exact_gaussian")
  expect_equal(test_object_params@gpfit_info$gpfit_parameters$gp_nugget, 1e-6)
  expect_true(length(test_object_params@gpfit_info$gp_symbolic_mean) > 0)
  expect_equal(test_object_params@gpfit_info$gp_symbolic_mean[[1]]$kernel, "gaussian")
  expect_true("alpha" %in% names(test_object_params@gpfit_info$gp_symbolic_mean[[1]]))
  expect_true("centers" %in% names(test_object_params@gpfit_info$gp_symbolic_mean[[1]]))
  expect_true("first_derivative" %in% names(test_object_params@gpfit_info$gp_symbolic_mean[[1]]))
  expect_true("second_derivative" %in% names(test_object_params@gpfit_info$gp_symbolic_mean[[1]]))
})

test_that("exact Gaussian GP analytic derivatives match finite differences", {
  timepoint <- seq(0, 1, length.out = 8)
  df <- data.frame(
    timepoint = timepoint,
    timepoint_n = seq_along(timepoint),
    od = sin(2 * pi * timepoint)
  )

  fit <- fit_predict_gp(
    df,
    t_steps = 8,
    gp_backend = "exact_gaussian",
    gp_lengthscale = 0.25,
    gp_nugget = 1e-6,
    optimize_hyperparameters = FALSE
  )

  eval_mean <- function(x) {
    pred <- .predict_exact_gaussian_gp(matrix(x, ncol = 1), fit)
    pred$prediction$mean
  }

  x <- matrix(seq(0.2, 0.8, length.out = 4), ncol = 1)
  h <- 1e-4
  analytic <- .exact_gaussian_derivatives(x, fit)
  numeric_first <- (eval_mean(x[, 1] + h) - eval_mean(x[, 1] - h)) / (2 * h)
  numeric_second <- (eval_mean(x[, 1] + h) - 2 * eval_mean(x[, 1]) + eval_mean(x[, 1] - h)) / h^2

  expect_equal(analytic$first_derivative, numeric_first, tolerance = 1e-5)
  expect_equal(analytic$second_derivative, numeric_second, tolerance = 1e-4)
})

test_that("exact Gaussian GP backend optimizes hyperparameters", {
  test_object_pre <- preprocess_data(test_object)

  fixed_fit <- fit_predict_gp(
    gather_od_data(test_object_pre, "Treatment", "Cancer_0.05"),
    t_steps = 10,
    gp_backend = "exact_gaussian",
    prepare_dataframe = TRUE,
    estimate_derivatives = TRUE,
    optimize_hyperparameters = FALSE
  )

  optimized_fit <- fit_predict_gp(
    gather_od_data(test_object_pre, "Treatment", "Cancer_0.05"),
    t_steps = 10,
    gp_backend = "exact_gaussian",
    prepare_dataframe = TRUE,
    estimate_derivatives = TRUE,
    optimize_hyperparameters = TRUE
  )

  expect_equal(optimized_fit$backend, "exact_gaussian")
  expect_false(isTRUE(all.equal(fixed_fit$lengthscale, optimized_fit$lengthscale)))
  expect_true(is.finite(optimized_fit$llik))
  expect_equal(optimized_fit$mle$conv, 0)
})

test_that("exact Gaussian GP defaults mirror laGP formulas without laGP calls", {
  X <- matrix(seq(0, 24, length.out = 145), ncol = 1)
  y <- sin(X[, 1])

  dist_sq <- as.matrix(stats::dist(X))^2
  positive_dist_sq <- dist_sq[upper.tri(dist_sq)]
  positive_dist_sq <- positive_dist_sq[positive_dist_sq > 0]
  r2s <- (y - mean(y))^2

  initial_params <- .initial_exact_gaussian_params(X, y)

  expect_equal(
    initial_params$lengthscale,
    as.numeric(stats::quantile(positive_dist_sq, probs = 0.1, names = FALSE))
  )
  expect_equal(initial_params$lengthscale_lower, min(positive_dist_sq) / 2)
  expect_equal(initial_params$lengthscale_upper, max(positive_dist_sq))
  expect_equal(
    initial_params$nugget_scale,
    max(
      as.numeric(stats::quantile(r2s, probs = 0.025, names = FALSE)),
      sqrt(.Machine$double.eps)
    )
  )
  expect_equal(initial_params$nugget_lower, sqrt(.Machine$double.eps))
  expect_equal(initial_params$nugget_upper, max(r2s))
})

test_that("old exact backend name maps to exact Gaussian", {
  test_object_pre <- preprocess_data(test_object)
  fit <- fit_predict_gp(
    gather_od_data(test_object_pre, "curve_id", metadata$curve_id[[1]]),
    t_steps = 10,
    gp_backend = "exact",
    prepare_dataframe = TRUE
  )

  expect_equal(fit$backend, "exact_gaussian")
})

test_that("subsample_timepoints keeps requested fraction with even spacing", {
  test_object_pre <- preprocess_data(test_object)

  thinned <- subsample_timepoints(
    test_object_pre,
    keep_fraction = 0.1,
    method = "even",
    always_keep = range(test_object_pre@od_data$timepoint)
  )

  expect_s4_class(thinned, "DGrowthR")
  expect_equal(
    length(unique(thinned@od_data$timepoint)),
    max(2, ceiling(length(unique(test_object_pre@od_data$timepoint)) * 0.1))
  )
  expect_true(all(range(test_object_pre@od_data$timepoint) %in% unique(thinned@od_data$timepoint)))
  expect_equal(nrow(thinned@metadata), nrow(test_object_pre@metadata))
})

test_that("subsample_timepoints defaults to GP information selection", {
  test_object_pre <- preprocess_data(test_object)

  thinned <- subsample_timepoints(
    test_object_pre,
    n_timepoints = 8,
    always_keep = range(test_object_pre@od_data$timepoint)
  )

  expect_s4_class(thinned, "DGrowthR")
  expect_equal(length(unique(thinned@od_data$timepoint)), 8)
  expect_true(all(range(test_object_pre@od_data$timepoint) %in% unique(thinned@od_data$timepoint)))
  expect_equal(attr(thinned@od_data, "subsampled_timepoints"), sort(unique(thinned@od_data$timepoint)))
})

test_that("subsample_timepoints can use growth-parameter landmarks", {
  test_object_pre <- preprocess_data(test_object)
  test_object_params <- estimate_growth_parameters(
    test_object_pre,
    model_covariate = "Treatment",
    predict_n_steps = 10,
    save_gp_data = FALSE
  )

  thinned <- subsample_timepoints(
    test_object_params,
    n_timepoints = 8,
    method = "growth_parameters",
    group_by = "Treatment",
    gpfit_id_column = "Treatment",
    always_keep = range(test_object_params@od_data$timepoint)
  )

  expect_s4_class(thinned, "DGrowthR")
  expect_equal(length(unique(thinned@od_data$timepoint)), 8)
  expect_equal(nrow(thinned@growth_parameters), 0)
  expect_true(all(range(test_object_params@od_data$timepoint) %in% unique(thinned@od_data$timepoint)))
})

test_that("subsample_timepoints can use derivative landmarks", {
  test_object_pre <- preprocess_data(test_object)
  test_object_params <- estimate_growth_parameters(
    test_object_pre,
    model_covariate = "Treatment",
    predict_n_steps = 10,
    save_gp_data = TRUE
  )

  thinned <- subsample_timepoints(
    test_object_params,
    n_timepoints = 8,
    method = "gp_derivative",
    group_by = "Treatment",
    gpfit_id_column = "Treatment",
    always_keep = range(test_object_params@od_data$timepoint)
  )

  expect_s4_class(thinned, "DGrowthR")
  expect_equal(length(unique(thinned@od_data$timepoint)), 8)
  expect_equal(length(thinned@gpfit_info), 0)
})

test_that("subsample_timepoints errors clearly when derivative data is missing", {
  test_object_pre <- preprocess_data(test_object)

  expect_error(
    subsample_timepoints(test_object_pre, method = "gp_derivative"),
    "requires object@gpfit_info\\$gpfit_data"
  )
})

laGP::deleteGPs()
laGP::deleteGPseps()
