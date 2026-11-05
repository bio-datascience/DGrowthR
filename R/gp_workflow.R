#' Fit Gaussian processes for the groups in a DGrowthR object
#'
#' @param object A DGrowthR object.
#' @param model_covariate Metadata column defining independent GP groups.
#' @param downsample_every_n_timepoints Use every nth observed timepoint.
#' @param gp_backend Either `"lagp"` or `"exact_gaussian"`.
#' @param gp_lengthscale Optional exact-GP lengthscale.
#' @param gp_nugget Optional exact-GP nugget.
#' @param gp_replicate_strategy Either `"raw"` or `"aggregate"`.
#' @param optimize_hyperparameters Optimize lengthscale and nugget.
#' @param n_cores Number of fitting workers.
#' @return A DGrowthR object with serializable fits in `gpfit_info$fits`.
#' @export
fit_gp <- function(object, model_covariate = "curve_id",
                   downsample_every_n_timepoints = 1,
                   gp_backend = c("lagp", "exact_gaussian", "exact"),
                   gp_lengthscale = NULL, gp_nugget = NULL,
                   gp_replicate_strategy = c("raw", "aggregate"),
                   optimize_hyperparameters = TRUE, n_cores = 1) {
  stopifnot(methods::is(object, "DGrowthR"))
  gp_backend <- .standardize_gp_backend(gp_backend)
  gp_replicate_strategy <- match.arg(gp_replicate_strategy)
  if (!model_covariate %in% names(object@metadata)) {
    stop("model_covariate is not present in metadata")
  }
  groups <- unique(object@metadata[[model_covariate]])
  if (anyNA(groups)) stop("model_covariate must not contain missing groups")
  fit_one <- function(group) {
    observations <- gather_od_data(object, model_covariate, group,
                                   downsample_every_n_timepoints)
    observations <- observations[, c("timepoint", "timepoint_n", "od"), drop = FALSE]
    if (nrow(observations) == 0L || .check_flat(observations)) {
      return(list(observations = observations, fit = NULL,
                  status = "flat_curve", error = NA_character_))
    }
    tryCatch({
      fit <- .fit_gp(observations, gp_backend = gp_backend,
                     gp_lengthscale = gp_lengthscale, gp_nugget = gp_nugget,
                     gp_replicate_strategy = gp_replicate_strategy,
                     optimize_hyperparameters = optimize_hyperparameters)
      if (identical(gp_backend, "lagp")) {
        laGP::deleteGPsep(fit$gpsepi)
        fit$gpsepi <- NULL
      }
      list(observations = observations, fit = fit, status = "fitted",
           error = NA_character_)
    }, error = function(e) {
      list(observations = observations, fit = NULL, status = "error",
           error = conditionMessage(e))
    })
  }
  if (n_cores > 1L) {
    cluster <- snow::makeSOCKcluster(n_cores)
    on.exit(snow::stopCluster(cluster), add = TRUE)
    doSNOW::registerDoSNOW(cluster)
    fits <- foreach::foreach(group = groups, .packages = "DGrowthR") %dopar% fit_one(group)
  } else {
    fits <- lapply(groups, fit_one)
  }
  names(fits) <- as.character(groups)
  object@gpfit_info <- list(
    fits = fits,
    gpfit_parameters = list(model_covariate = model_covariate,
                            subsampled_timepoints = attr(object@od_data, "subsampled_timepoints"),
                            downsample_every_n_timepoints = downsample_every_n_timepoints,
                            gp_backend = gp_backend, gp_lengthscale = gp_lengthscale,
                            gp_nugget = gp_nugget,
                            gp_replicate_strategy = gp_replicate_strategy,
                            optimize_hyperparameters = optimize_hyperparameters),
    gp_symbolic_mean = lapply(fits, function(entry) {
      if (is.null(entry$fit)) NULL else entry$fit$symbolic_mean
    }))
  object@growth_parameters <- data.frame()
  object
}

.predict_saved_gp <- function(fit, timepoints, complete_sigma = FALSE,
                              estimate_derivatives = TRUE) {
  if (identical(fit$backend, "lagp")) {
    fit$gpsepi <- laGP::newGPsep(fit$X, fit$y, d = fit$lengthscale,
                                 g = fit$nugget, dK = FALSE)
  }
  .predict_gp(matrix(timepoints, ncol = 1,
                     dimnames = list(NULL, "timepoint")), fit,
              complete_sigma = complete_sigma, prepare_dataframe = TRUE,
              estimate_derivatives = estimate_derivatives, delete = TRUE)
}

#' Predict from stored Gaussian-process fits
#'
#' @param object A DGrowthR object returned by [fit_gp()].
#' @param predict_n_steps Number of equally spaced prediction timepoints.
#' @param prediction_timepoints Optional explicit numeric prediction timepoints.
#' @param complete_sigma Store the full latent posterior covariance.
#' @param estimate_derivatives Calculate first and second derivatives.
#' @return A DGrowthR object containing `gpfit_info$gpfit_data`.
#' @export
predict_gp <- function(object, predict_n_steps = 100,
                       prediction_timepoints = NULL, complete_sigma = FALSE,
                       estimate_derivatives = TRUE) {
  stopifnot(methods::is(object, "DGrowthR"))
  if (is.null(object@gpfit_info$fits)) stop("Call fit_gp() first")
  if (!is.null(prediction_timepoints) &&
      (!is.numeric(prediction_timepoints) || anyNA(prediction_timepoints) ||
       any(!is.finite(prediction_timepoints)))) {
    stop("prediction_timepoints must be finite numeric values")
  }
  if (is.null(prediction_timepoints) &&
      (length(predict_n_steps) != 1L || predict_n_steps < 2L)) {
    stop("predict_n_steps must be at least 2")
  }
  frames <- list()
  for (id in names(object@gpfit_info$fits)) {
    entry <- object@gpfit_info$fits[[id]]
    if (is.null(entry$fit)) next
    times <- if (is.null(prediction_timepoints)) {
      seq(min(entry$observations$timepoint), max(entry$observations$timepoint),
          length.out = predict_n_steps)
    } else sort(unique(prediction_timepoints))
    predicted <- .predict_saved_gp(entry$fit, times, complete_sigma,
                                   estimate_derivatives)
    entry$prediction <- predicted$prediction
    object@gpfit_info$fits[[id]] <- entry
    frame <- predicted$prediction$prediction_dataframe
    frame$gpfit_id <- id
    frames[[id]] <- frame
  }
  object@gpfit_info$gpfit_data <- dplyr::bind_rows(frames)
  object@gpfit_info$gpfit_parameters$predict_n_steps <- predict_n_steps
  object@gpfit_info$gpfit_parameters$prediction_timepoints <- prediction_timepoints
  object@gpfit_info$gpfit_parameters$complete_sigma <- complete_sigma
  object
}

#' Assess fit quality at observed timepoints
#'
#' @param object A DGrowthR object returned by [fit_gp()].
#' @param level Two-sided interval probability.
#' @param interval_type Use `"predictive"` for group-mean observations or
#'   `"latent"` for the latent GP mean.
#' @param review_duration,review_points Thresholds for a review flag.
#' @param fail_duration,fail_excess Joint thresholds for a fail flag.
#' @return A list with `summary`, `runs`, and `points` data frames.
#' @export
assess_gp_fit <- function(object, level = 0.95,
                          interval_type = c("predictive", "latent"),
                          review_duration = 0.5, review_points = 4L,
                          fail_duration = 1, fail_excess = 1) {
  stopifnot(methods::is(object, "DGrowthR"))
  if (is.null(object@gpfit_info$fits)) stop("Call fit_gp() first")
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must lie strictly between zero and one")
  }
  interval_type <- match.arg(interval_type)
  summaries <- runs <- points <- list()
  for (id in names(object@gpfit_info$fits)) {
    entry <- object@gpfit_info$fits[[id]]
    if (is.null(entry$fit)) {
      summaries[[id]] <- data.frame(gpfit_id = id, status = entry$status,
                                    optimizer_convergence = NA_integer_)
      next
    }
    obs <- entry$observations %>%
      group_by(timepoint) %>%
      summarise(observed = mean(od, na.rm = TRUE),
                n_replicates = sum(is.finite(od)), .groups = "drop") %>%
      arrange(timepoint)
    times <- obs$timepoint
    fit <- entry$fit
    if (identical(fit$backend, "lagp")) {
      gp <- laGP::newGPsep(fit$X, fit$y, d = fit$lengthscale,
                           g = fit$nugget, dK = FALSE)
      latent <- laGP::predGPsep(gp, matrix(times, ncol = 1),
                                lite = TRUE, nonug = TRUE)
      noisy <- laGP::predGPsep(gp, matrix(times, ncol = 1),
                               lite = TRUE, nonug = FALSE)
      laGP::deleteGPsep(gp)
      quantile <- stats::qt((1 + level) / 2, df = latent$df)
      latent_sd <- sqrt(pmax(latent$s2, 0))
      noise_variance <- pmax(noisy$s2 - latent$s2, 0)
      predictive_sd <- sqrt(pmax(latent$s2 + noise_variance /
                                   pmax(obs$n_replicates, 1), 0))
      mu <- latent$mean
    } else {
      predicted <- .predict_exact_gaussian_gp(matrix(times, ncol = 1), fit,
                                               complete_sigma = FALSE)$prediction
      mu <- predicted$mean
      quantile <- stats::qnorm((1 + level) / 2)
      latent_sd <- sqrt(pmax(predicted$s2, 0))
      predictive_sd <- sqrt(pmax(predicted$s2 + fit$sigma2 * fit$nugget /
                                   pmax(obs$n_replicates, 1), 0))
    }
    obs$mean <- as.numeric(mu)
    obs$latent_lower <- mu - quantile * latent_sd
    obs$latent_upper <- mu + quantile * latent_sd
    obs$predictive_lower <- mu - quantile * predictive_sd
    obs$predictive_upper <- mu + quantile * predictive_sd
    lower <- obs[[paste0(interval_type, "_lower")]]
    upper <- obs[[paste0(interval_type, "_upper")]]
    direction <- ifelse(obs$observed > upper, "above",
                        ifelse(obs$observed < lower, "below", "inside"))
    width <- pmax(upper - lower, sqrt(.Machine$double.eps))
    excess <- ifelse(direction == "above", (obs$observed - upper) / width,
                     ifelse(direction == "below", (lower - obs$observed) / width, 0))
    obs$direction <- direction
    obs$relative_excess <- excess
    obs$gpfit_id <- id
    obs$interval_type <- interval_type
    points[[id]] <- obs
    breaks <- cumsum(c(TRUE, direction[-1L] != direction[-length(direction)]))
    run_table <- dplyr::bind_rows(lapply(split(seq_along(times), breaks), function(ix) {
      if (direction[ix[1L]] == "inside") return(NULL)
      data.frame(gpfit_id = id, direction = direction[ix[1L]],
                 start = times[ix[1L]], end = times[tail(ix, 1L)],
                 duration = times[tail(ix, 1L)] - times[ix[1L]],
                 n_points = length(ix), max_relative_excess = max(excess[ix]),
                 mean_relative_excess = mean(excess[ix]))
    }))
    runs[[id]] <- run_table
    residual <- obs$observed - mu
    observed_auc <- pracma::trapz(times, obs$observed)
    fitted_auc <- pracma::trapz(times, mu)
    longest <- if (nrow(run_table)) run_table[which.max(run_table$duration), ] else NULL
    flag <- if (nrow(run_table) && any(run_table$duration >= fail_duration &
                                      run_table$max_relative_excess >= fail_excess)) {
      "fail"
    } else if (nrow(run_table) && any(run_table$duration >= review_duration |
                                       run_table$n_points >= review_points)) {
      "review"
    } else "pass"
    if (!is.na(fit$mle$conv) && fit$mle$conv != 0L && flag == "pass") {
      flag <- "review"
    }
    summaries[[id]] <- data.frame(
      gpfit_id = id, status = flag,
      optimizer_convergence = fit$mle$conv,
      rmse = sqrt(mean(residual^2)), mae = mean(abs(residual)),
      normalized_rmse = sqrt(mean(residual^2)) /
        pmax(diff(range(obs$observed)), sqrt(.Machine$double.eps)),
      coverage = mean(direction == "inside"), outside_fraction = mean(direction != "inside"),
      observed_auc = observed_auc, fitted_auc = fitted_auc,
      relative_auc_error = abs(observed_auc - fitted_auc) /
        pmax(abs(observed_auc), sqrt(.Machine$double.eps)),
      longest_outside_duration = if (is.null(longest)) 0 else longest$duration,
      longest_outside_points = if (is.null(longest)) 0L else longest$n_points,
      longest_outside_start = if (is.null(longest)) NA_real_ else longest$start,
      longest_outside_end = if (is.null(longest)) NA_real_ else longest$end,
      longest_outside_direction = if (is.null(longest)) NA_character_ else longest$direction,
      max_relative_excess = if (nrow(run_table)) max(run_table$max_relative_excess) else 0,
      mean_relative_excess = if (any(excess > 0)) mean(excess[excess > 0]) else 0,
      total_outside_duration = if (nrow(run_table)) sum(run_table$duration) else 0,
      interval_type = interval_type, interval_level = level)
  }
  list(summary = dplyr::bind_rows(summaries), runs = dplyr::bind_rows(runs),
       points = dplyr::bind_rows(points))
}
