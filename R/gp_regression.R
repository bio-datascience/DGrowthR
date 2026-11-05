#------------------------------------------------------------------------------------------------------------------
#' Fit a Gaussian process model to data with laGP
#'
#' This function fits a Gaussian process model to the provided data. It is an internal function that
#' really only does the fit. Nothing else is done or returned by this function.
#'
#' @param X A data.frame containing the covariates used for model fitting.
#' @param y A vector containing the values to be predicted.
#' @param delete A logical value. If TRUE, the Gaussian process model will be deleted after fitting.
#' @param optimize_hyperparameters A logical value. If TRUE, estimate the
#' lengthscale and nugget by maximum likelihood. If FALSE, use the initial
#' values from `laGP::darg()` and `laGP::garg()` without optimization.
#'
#' @return The input object, with the Gaussian process information as a list.
#'
#' @keywords internal
.fit_lagp <- function(X, y, delete=FALSE, optimize_hyperparameters=TRUE) {

  lagp_defaults <- .lagp_hyperparameter_defaults(X, y)
  d <- lagp_defaults$d
  g <- lagp_defaults$g

  #print(g)
  # Start an anisotropic GP model
  gpsepi <- laGP::newGPsep(X, y, d = d$start, g = g$start, dK = TRUE)
  # Fit GP using MLE to estimate the best lengthscale and nugget
  if(optimize_hyperparameters){
    mle <- laGP::mleGPsep(gpsepi, param = "both", tmin = c(d$min, g$min), tmax = c(d$max, g$max))
  }else{
    mle <- list(conv = NA_integer_,
                msg = "Hyperparameter optimization skipped")
  }
  # Compute the log-likelihood
  llik <- laGP::llikGPsep(gpsepi)
  # Compute posterior
  posterior <- laGP::llikGPsep(gpsepi, dab = d$ab, gab = g$ab)

  # Do not include gpsep obj by default
  theta <- if (optimize_hyperparameters) mle$theta else c(d$start, g$start)
  fit_info <- list(X = X, y = y, mle = mle, llik = llik,
                   posterior = posterior, backend = "lagp",
                   lengthscale = theta[seq_len(ncol(X))],
                   nugget = tail(theta, 1))

  if (delete == TRUE) {
    laGP::deleteGPsep(gpsepi)
    fit_info$gpsepi <- NULL
  } else if (delete == FALSE) {
    fit_info$gpsepi <- gpsepi
  }

  return(fit_info)
}


#------------------------------------------------------------------------------------------------------------------
#' Create laGP hyperparameter defaults
#'
#' @param X A matrix containing the covariates used for model fitting.
#' @param y A vector containing the values to be predicted.
#'
#' @return A list with laGP-style range (`d`) and nugget (`g`) arguments.
#'
#' @keywords internal
.lagp_hyperparameter_defaults <- function(X, y) {

  d <- laGP::darg(NULL, X)

  g <- tryCatch({
    laGP::garg(list(mle = TRUE), y)
  },
  error=function(cond){
    message("DGrowthR is setting new nugget priors...")

    y_var <- stats::var(y)
    if(is.na(y_var) || y_var <= 0){
      y_var <- sqrt(.Machine$double.eps)
    }

    r2mean <- mean((y - mean(y))^2)

    if(r2mean == 0){
      r2mean <- median((y - mean(y))^2)
    }

    if(r2mean == 0){
      r2mean <- max((y - mean(y))^2)
    }

    if(r2mean == 0 || is.na(r2mean)){
      r2mean <- sqrt(.Machine$double.eps)
    }

    list(mle=TRUE,
         start = max(0.1 * y_var, sqrt(.Machine$double.eps)),
         min = sqrt(.Machine$double.eps),
         max = y_var,
         ab=c(3/2, 1276.323 / r2mean))
  })

  list(d = d, g = g)
}


#------------------------------------------------------------------------------------------------------------------
#' Fit a lightweight exact Gaussian process model
#'
#' @param X A matrix containing the covariates used for model fitting.
#' @param y A vector containing the values to be predicted.
#' @param lengthscale Optional Gaussian range value or vector, one value per
#' covariate. If NULL, a robust default is estimated from the input range.
#' @param nugget Optional relative observation-noise variance used as a nugget.
#' If NULL, a robust default is used.
#' @param optimize_hyperparameters A logical value. If TRUE, estimate the
#' Gaussian range, process variance, and nugget by maximum marginal likelihood.
#'
#' @return A list with the fitted exact GP state.
#'
#' @keywords internal
.fit_exact_gaussian_gp <- function(X, y, lengthscale=NULL, nugget=NULL,
                                   optimize_hyperparameters=TRUE) {

  X <- as.matrix(X)
  dist_sq <- .gaussian_squared_distances(X, X)
  y_mean <- 0
  y_centered <- y - y_mean
  initial_params <- .initial_exact_gaussian_params(X, y, lengthscale, nugget)

  if(optimize_hyperparameters){
    opt_start <- .exact_gaussian_optimization_starts(initial_params)
    opt_lower <- log(c(initial_params$lengthscale_lower,
                       initial_params$nugget_lower))
    opt_upper <- log(c(initial_params$lengthscale_upper,
                       initial_params$nugget_upper))
    opt <- .exact_gaussian_profile_optimize(opt_start, opt_lower, opt_upper,
                                            X = X,
                                            dist_sq = dist_sq,
                                            y_centered = y_centered,
                                            nugget_weights = initial_params$nugget_weights)

    fit_params <- .decode_exact_gaussian_profile_params(opt$par, ncol(X))
    sigma2_fit_state <- .exact_gaussian_fit_state(X,
                                                  y_centered,
                                                  fit_params$lengthscale,
                                                  sigma2 = 1,
                                                  nugget = fit_params$nugget * initial_params$nugget_weights,
                                                  dist_sq = dist_sq)
    fit_params$sigma2 <- .profile_exact_gaussian_sigma2(y_centered,
                                                        sigma2_fit_state$alpha)
    opt_message <- opt$message
    if(is.null(opt_message)){
      opt_message <- "Exact Gaussian GP lengthscale and nugget optimized; sigma2 profiled"
    }

    mle <- list(conv = opt$convergence,
                msg = opt_message,
                value = opt$value,
                counts = opt$counts)
  }else{
    fit_params <- list(lengthscale = initial_params$lengthscale,
                       sigma2 = initial_params$sigma2,
                       nugget = initial_params$nugget_scale)
    mle <- list(conv = NA_integer_,
                msg = "Exact Gaussian GP uses fixed hyperparameters")
  }

  fit_state <- .exact_gaussian_fit_state(X,
                                         y_centered,
                                         fit_params$lengthscale,
                                         fit_params$sigma2,
                                         fit_params$nugget * initial_params$nugget_weights,
                                         dist_sq = dist_sq)

  llik <- -fit_state$neg_loglik

  list(X = X,
       y = y,
       y_mean = y_mean,
       sigma2 = fit_params$sigma2,
       lengthscale = fit_params$lengthscale,
       nugget = fit_params$nugget,
       nugget_weights = initial_params$nugget_weights,
       jitter = fit_state$jitter,
       chol = fit_state$chol,
       alpha = fit_state$alpha,
       symbolic_mean = .exact_gaussian_symbolic_mean(X, y_mean, fit_params$sigma2,
                                                     fit_params$lengthscale,
                                                     fit_state$alpha),
       mle = mle,
       llik = llik,
       posterior = llik,
       backend = "exact_gaussian")
}


#------------------------------------------------------------------------------------------------------------------
#' Optimize exact Gaussian GP profiled likelihood
#'
#' @param start Initial log-parameters.
#' @param lower Lower log-parameter bounds.
#' @param upper Upper log-parameter bounds.
#' @param X A matrix containing the covariates used for model fitting.
#' @param dist_sq Precomputed per-covariate squared distances.
#' @param y_centered A centered response vector.
#' @param nugget_weights Optional relative nugget weights.
#'
#' @return An optimizer result list.
#'
#' @keywords internal
.exact_gaussian_profile_optimize <- function(start, lower, upper, X, dist_sq,
                                             y_centered, nugget_weights) {

  if(is.null(dim(start))){
    starts <- matrix(start, nrow = 1)
  }else{
    starts <- start
  }

  best_opt <- NULL
  for(i in seq_len(nrow(starts))){
    current_start <- starts[i, ]

    opt <- NULL
    if(is.loaded("DGrowthR_exact_gaussian_lbfgsb_dist")){
      opt <- .Call("DGrowthR_exact_gaussian_lbfgsb_dist",
                   as.numeric(current_start),
                   as.numeric(lower),
                   as.numeric(upper),
                   dist_sq,
                   as.numeric(y_centered),
                   as.numeric(nugget_weights),
                   as.integer(100))
      opt$message <- paste("C L-BFGS-B:", opt$message)
    }

    if(is.null(opt) || !identical(opt$convergence, 0L)){
      opt <- stats::optim(
        par = current_start,
        fn = .exact_gaussian_profile_neg_loglik,
        X = X,
        dist_sq = dist_sq,
        y_centered = y_centered,
        nugget_weights = nugget_weights,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper
      )
    }

    if(is.null(best_opt) || opt$value < best_opt$value){
      best_opt <- opt
      best_opt$start_index <- i
    }
  }

  if(!is.null(best_opt$message)){
    best_opt$message <- paste0(best_opt$message, "; selected multistart ", best_opt$start_index)
  }else{
    best_opt$message <- paste("selected multistart", best_opt$start_index)
  }

  best_opt
}


#------------------------------------------------------------------------------------------------------------------
#' Create deterministic optimization starts for an exact Gaussian GP
#'
#' @param initial_params Output from `.initial_exact_gaussian_params()`.
#'
#' @return A matrix of log-scale optimizer starts.
#'
#' @keywords internal
.exact_gaussian_optimization_starts <- function(initial_params) {

  lower <- c(initial_params$lengthscale_lower, initial_params$nugget_lower)
  start <- c(initial_params$lengthscale, initial_params$nugget_scale)
  upper <- c(initial_params$lengthscale_upper, initial_params$nugget_upper)

  small_lengthscale <- initial_params$lengthscale_lower * 8
  mid_lengthscale <- sqrt(initial_params$lengthscale_lower * initial_params$lengthscale_upper)
  large_lengthscale <- sqrt(initial_params$lengthscale * initial_params$lengthscale_upper)
  low_nugget <- pmax(initial_params$nugget_lower * 10, sqrt(.Machine$double.eps))
  mid_nugget <- sqrt(initial_params$nugget_lower * initial_params$nugget_upper)
  high_nugget <- sqrt(initial_params$nugget_scale * initial_params$nugget_upper)

  starts <- rbind(
    start,
    c(small_lengthscale, initial_params$nugget_scale),
    c(mid_lengthscale, initial_params$nugget_scale),
    c(large_lengthscale, initial_params$nugget_scale),
    c(initial_params$lengthscale, low_nugget),
    c(initial_params$lengthscale, mid_nugget),
    c(initial_params$lengthscale, high_nugget)
  )

  starts <- sweep(starts, 2, lower, pmax)
  starts <- sweep(starts, 2, upper, pmin)
  starts <- unique(signif(starts, digits = 12))
  log(starts)
}


#------------------------------------------------------------------------------------------------------------------
#' Create initial hyperparameters for an exact Gaussian GP
#'
#' @param X A matrix containing the covariates used for model fitting.
#' @param y A vector containing the values to be predicted.
#' @param lengthscale Optional user-specified lengthscale.
#' @param nugget Optional user-specified nugget.
#'
#' @return A list with initial values and optimizer upper bounds.
#'
#' @keywords internal
.initial_exact_gaussian_params <- function(X, y, lengthscale=NULL, nugget=NULL) {

  X <- as.matrix(X)
  user_lengthscale <- !is.null(lengthscale)
  user_nugget <- !is.null(nugget)
  lengthscale_defaults <- .exact_gaussian_lengthscale_defaults(X)
  y_var <- stats::var(y)

  if(is.na(y_var) || y_var <= 0){
    y_var <- sqrt(.Machine$double.eps)
  }

  if(is.null(nugget)){
    r2s <- (y - mean(y))^2
    nugget_scale <- as.numeric(stats::quantile(r2s, probs = 0.025, names = FALSE))
    nugget_scale <- max(nugget_scale, sqrt(.Machine$double.eps))
    nugget_lower <- sqrt(.Machine$double.eps)
    nugget_upper <- max(r2s, na.rm = TRUE)
    nugget_upper <- max(nugget_upper, nugget_scale, 10 * nugget_lower)
    nugget_weights <- rep(1, nrow(X))
  }else if(length(nugget) == 1){
    nugget_scale <- nugget
    nugget_lower <- min(sqrt(.Machine$double.eps), nugget)
    nugget_upper <- max(nugget_scale, 10 * nugget)
    nugget_weights <- rep(1, nrow(X))
  }else if(length(nugget) == nrow(X)){
    nugget_scale <- 1
    nugget_lower <- sqrt(.Machine$double.eps)
    nugget_upper <- 100
    nugget_weights <- nugget
  }else{
    stop("nugget must be NULL, length 1, or length nrow(X).")
  }

  if(is.null(lengthscale)){
    lengthscale <- lengthscale_defaults$start
  }

  if(length(lengthscale) == 1){
    lengthscale <- rep(lengthscale, ncol(X))
  }

  lengthscale_lower <- lengthscale_defaults$lower
  lengthscale_upper <- lengthscale_defaults$upper
  if(user_lengthscale){
    lengthscale_lower <- pmin(lengthscale_lower, lengthscale)
    lengthscale_upper <- pmax(lengthscale_upper, 10 * lengthscale)
  }

  list(lengthscale = lengthscale,
       sigma2 = y_var,
       nugget_scale = nugget_scale,
       nugget_weights = nugget_weights,
       lengthscale_lower = lengthscale_lower,
       lengthscale_upper = lengthscale_upper,
       sigma2_upper = max(100 * y_var, 1),
       nugget_lower = nugget_lower,
       nugget_upper = nugget_upper)
}


#------------------------------------------------------------------------------------------------------------------
#' Create exact Gaussian GP range defaults without laGP
#'
#' @param X A matrix containing the covariates used for model fitting.
#'
#' @return A list with start, lower, and upper range values.
#'
#' @keywords internal
.exact_gaussian_lengthscale_defaults <- function(X) {

  X <- as.matrix(X)
  if(nrow(X) < 2){
    return(list(start = rep(1, ncol(X)),
                lower = rep(sqrt(.Machine$double.eps), ncol(X)),
                upper = rep(100, ncol(X))))
  }

  dist_sq <- .gaussian_squared_distances(X, X)
  if(length(dim(dist_sq)) == 3){
    dist_sq <- apply(dist_sq, c(1, 2), sum)
  }

  distances <- dist_sq[upper.tri(dist_sq)]
  distances <- distances[is.finite(distances) & distances > 0]

  if(length(distances) == 0){
    return(list(start = rep(1, ncol(X)),
                lower = rep(sqrt(.Machine$double.eps), ncol(X)),
                upper = rep(100, ncol(X))))
  }

  d_start <- as.numeric(stats::quantile(distances, probs = 0.1, names = FALSE))
  d_min <- min(distances, na.rm = TRUE) / 2
  d_max <- max(distances, na.rm = TRUE)

  d_start <- max(d_start, sqrt(.Machine$double.eps))
  d_min <- max(d_min, sqrt(.Machine$double.eps))
  d_max <- max(d_max, d_start, 10 * d_min)

  list(start = rep(d_start, ncol(X)),
       lower = rep(d_min, ncol(X)),
       upper = rep(d_max, ncol(X)))
}


#------------------------------------------------------------------------------------------------------------------
#' Decode exact Gaussian GP parameters from log scale
#'
#' @param log_params Parameters on the log scale.
#' @param n_covariates Number of covariates.
#'
#' @return A list with lengthscale, sigma2, and nugget.
#'
#' @keywords internal
.decode_exact_gaussian_params <- function(log_params, n_covariates) {

  params <- exp(log_params)
# The "lengthscale" field follows laGP's Gaussian range parameterization:
# k(x, x') = sigma2 * exp(-sum((x - x')^2 / lengthscale)).
  list(lengthscale = params[seq_len(n_covariates)],
       sigma2 = params[n_covariates + 1],
       nugget = params[n_covariates + 2])
}


#------------------------------------------------------------------------------------------------------------------
#' Decode exact Gaussian GP profiled parameters from log scale
#'
#' @param log_params Parameters on the log scale.
#' @param n_covariates Number of covariates.
#'
#' @return A list with lengthscale and nugget.
#'
#' @keywords internal
.decode_exact_gaussian_profile_params <- function(log_params, n_covariates) {

  params <- exp(log_params)
  list(lengthscale = params[seq_len(n_covariates)],
       nugget = params[n_covariates + 1])
}


#------------------------------------------------------------------------------------------------------------------
#' Calculate exact Gaussian GP fit state and negative log-likelihood
#'
#' @param X A matrix containing the covariates used for model fitting.
#' @param y_centered A centered response vector.
#' @param lengthscale Gaussian-kernel lengthscale value or vector.
#' @param sigma2 Process variance.
#' @param nugget Relative observation-noise variance.
#'
#' @return A list with Cholesky factor, alpha, jitter, and negative
#' log-likelihood.
#'
#' @keywords internal
.exact_gaussian_fit_state <- function(X, y_centered, lengthscale, sigma2, nugget,
                                      dist_sq=NULL) {

  if(!is.null(dist_sq) && is.loaded("DGrowthR_exact_gaussian_fit_state_dist")){
    return(.Call("DGrowthR_exact_gaussian_fit_state_dist",
                 dist_sq,
                 as.numeric(y_centered),
                 as.numeric(lengthscale),
                 as.numeric(sigma2),
                 as.numeric(nugget)))
  }

  if(is.loaded("DGrowthR_exact_gaussian_fit_state")){
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    return(.Call("DGrowthR_exact_gaussian_fit_state",
                 X,
                 as.numeric(y_centered),
                 as.numeric(lengthscale),
                 as.numeric(sigma2),
                 as.numeric(nugget)))
  }

  .exact_gaussian_fit_state_r(X, y_centered, lengthscale, sigma2, nugget,
                              dist_sq = dist_sq)
}


#------------------------------------------------------------------------------------------------------------------
#' R fallback for exact Gaussian GP fit state
#'
#' @param X A matrix containing the covariates used for model fitting.
#' @param y_centered A centered response vector.
#' @param lengthscale Gaussian-kernel lengthscale value or vector.
#' @param sigma2 Process variance.
#' @param nugget Relative observation-noise variance.
#'
#' @return A list with Cholesky factor, alpha, jitter, and negative
#' log-likelihood.
#'
#' @keywords internal
.exact_gaussian_fit_state_r <- function(X, y_centered, lengthscale, sigma2, nugget,
                                        dist_sq=NULL) {

  K <- sigma2 * .gaussian_kernel_from_distances(X, X, lengthscale,
                                                dist_sq = dist_sq)
  diag(K) <- diag(K) + sigma2 * nugget
  chol_info <- .chol_with_jitter(K)
  alpha <- backsolve(chol_info$chol, forwardsolve(t(chol_info$chol), y_centered))
  neg_loglik <- 0.5 * sum(y_centered * alpha) +
    sum(log(diag(chol_info$chol))) +
    0.5 * length(y_centered) * log(2 * pi)

  list(chol = chol_info$chol,
       alpha = alpha,
       jitter = chol_info$jitter,
       neg_loglik = neg_loglik)
}


#------------------------------------------------------------------------------------------------------------------
#' Negative log-likelihood for exact Gaussian GP optimization
#'
#' @param log_params Parameters on the log scale.
#' @param X A matrix containing the covariates used for model fitting.
#' @param y_centered A centered response vector.
#'
#' @return Numeric negative log-likelihood.
#'
#' @keywords internal
.exact_gaussian_neg_loglik <- function(log_params, X, y_centered, nugget_weights=NULL,
                                       dist_sq=NULL) {

  params <- .decode_exact_gaussian_params(log_params, ncol(X))
  if(is.null(nugget_weights)){
    nugget_weights <- rep(1, nrow(X))
  }
  out <- tryCatch(
    .exact_gaussian_fit_state(X, y_centered, params$lengthscale,
                              params$sigma2, params$nugget * nugget_weights,
                              dist_sq = dist_sq)$neg_loglik,
    error = function(cond) Inf
  )

  if(!is.finite(out)) Inf else out
}


#------------------------------------------------------------------------------------------------------------------
#' Negative log-likelihood with profiled process variance
#'
#' @param log_params Lengthscale(s) and nugget on the log scale.
#' @param X A matrix containing the covariates used for model fitting.
#' @param y_centered A centered response vector.
#' @param nugget_weights Optional relative nugget weights.
#' @param dist_sq Optional precomputed per-covariate squared distances.
#'
#' @return Numeric negative log-likelihood with sigma2 profiled out.
#'
#' @keywords internal
.exact_gaussian_profile_neg_loglik <- function(log_params, X, y_centered,
                                               nugget_weights=NULL,
                                               dist_sq=NULL) {

  params <- .decode_exact_gaussian_profile_params(log_params, ncol(X))
  if(is.null(nugget_weights)){
    nugget_weights <- rep(1, nrow(X))
  }

  out <- tryCatch({
    fit_state <- .exact_gaussian_fit_state(X, y_centered, params$lengthscale,
                                           sigma2 = 1,
                                           nugget = params$nugget * nugget_weights,
                                           dist_sq = dist_sq)
    sigma2_hat <- .profile_exact_gaussian_sigma2(y_centered, fit_state$alpha)
    logdet_half <- sum(log(diag(fit_state$chol)))
    0.5 * length(y_centered) * (log(sigma2_hat) + 1 + log(2 * pi)) +
      logdet_half
  }, error = function(cond) Inf)

  if(!is.finite(out)) Inf else out
}


#------------------------------------------------------------------------------------------------------------------
#' Profile the exact Gaussian GP process variance
#'
#' @param y_centered A centered response vector.
#' @param correlation_alpha Inverse correlation matrix times y_centered.
#'
#' @return The profiled process variance.
#'
#' @keywords internal
.profile_exact_gaussian_sigma2 <- function(y_centered, correlation_alpha) {

  sigma2 <- sum(y_centered * correlation_alpha) / length(y_centered)
  max(sigma2, sqrt(.Machine$double.eps))
}


#------------------------------------------------------------------------------------------------------------------
#' Predict from a lightweight exact Gaussian process model
#'
#' @param XX A matrix containing prediction covariates.
#' @param exact_fit The output of `.fit_exact_gp()`.
#' @param complete_sigma A logical value indicating whether the full covariance
#' matrix should be returned.
#'
#' @return The input fit with a prediction slot.
#'
#' @keywords internal
.predict_exact_gaussian_gp <- function(XX, exact_fit, complete_sigma=FALSE) {

  XX <- as.matrix(XX)
  K_x_xnew <- exact_fit$sigma2 *
    .gaussian_range_kernel(exact_fit$X, XX, exact_fit$lengthscale)
  mean_pred <- as.vector(exact_fit$y_mean + t(K_x_xnew) %*% exact_fit$alpha)
  v <- forwardsolve(t(exact_fit$chol), K_x_xnew)

  if(complete_sigma){
    K_new_new <- exact_fit$sigma2 *
      .gaussian_range_kernel(XX, XX, exact_fit$lengthscale)
    Sigma <- K_new_new - crossprod(v)
    Sigma <- (Sigma + t(Sigma)) / 2
    diag(Sigma) <- pmax(diag(Sigma), sqrt(.Machine$double.eps))
    prediction <- list(mean = mean_pred,
                       Sigma = Sigma,
                       df = length(exact_fit$y))
  }else{
    s2 <- pmax(exact_fit$sigma2 - colSums(v^2), sqrt(.Machine$double.eps))
    prediction <- list(mean = mean_pred,
                       s2 = s2,
                       df = length(exact_fit$y))
  }

  exact_fit$prediction <- prediction
  exact_fit
}


#------------------------------------------------------------------------------------------------------------------
#' Analytic derivatives from an exact Gaussian-kernel posterior mean
#'
#' @param XX A matrix containing prediction covariates.
#' @param exact_fit The output of `.fit_exact_gp()`.
#' @param variable_index Index of the covariate to differentiate with respect
#' to. DGrowthR uses timepoint as the first modelling covariate.
#'
#' @return A data.frame with first_derivative and second_derivative columns.
#'
#' @keywords internal
.exact_gaussian_derivatives <- function(XX, exact_fit, variable_index=1) {

  XX <- as.matrix(XX)
  K_x_xnew <- exact_fit$sigma2 *
    .gaussian_range_kernel(exact_fit$X, XX, exact_fit$lengthscale)
  diff_mat <- outer(exact_fit$X[, variable_index], XX[, variable_index], "-")
  lengthscale <- exact_fit$lengthscale[variable_index]
  weighted_kernel <- sweep(K_x_xnew, 1, exact_fit$alpha, "*")

  first_derivative <- colSums(weighted_kernel * (2 * diff_mat / lengthscale))
  second_derivative <- colSums(
    weighted_kernel * ((4 * diff_mat^2 / lengthscale^2) - (2 / lengthscale))
  )

  data.frame(first_derivative = first_derivative,
             second_derivative = second_derivative)
}


#------------------------------------------------------------------------------------------------------------------
#' Symbolic posterior mean for an exact Gaussian-kernel GP
#'
#' @param X A matrix containing the fitted covariates.
#' @param y_mean The fitted constant mean.
#' @param sigma2 The fitted process variance.
#' @param lengthscale Gaussian-kernel lengthscale value or vector.
#' @param alpha The fitted GP weights.
#'
#' @return A list with a symbolic formula template and fitted constants.
#'
#' @keywords internal
.exact_gaussian_symbolic_mean <- function(X, y_mean, sigma2, lengthscale, alpha) {

  variable_names <- colnames(X)
  if(is.null(variable_names)){
    variable_names <- paste0("x", seq_len(ncol(X)))
  }

  list(
    kernel = "gaussian",
    formula = "f(x) = y_mean + sum_i alpha_i * sigma2 * exp(-sum_j ((x_j - X_ij)^2 / lengthscale_j))",
    first_derivative = "df/dx_k = sum_i alpha_i * k(x, X_i) * 2 * (X_ik - x_k) / lengthscale_k",
    second_derivative = "d2f/dx_k2 = sum_i alpha_i * k(x, X_i) * ((4 * (X_ik - x_k)^2 / lengthscale_k^2) - (2 / lengthscale_k))",
    variable_names = variable_names,
    centers = as.data.frame(X),
    alpha = alpha,
    y_mean = y_mean,
    sigma2 = sigma2,
    lengthscale = lengthscale
  )
}


#------------------------------------------------------------------------------------------------------------------
#' Gaussian covariance kernel with laGP-style range parameterization
#'
#' @param X1 A numeric matrix.
#' @param X2 A numeric matrix.
#' @param lengthscale A numeric vector of Gaussian range parameters.
#'
#' @return A covariance matrix.
#'
#' @keywords internal
.gaussian_range_kernel <- function(X1, X2, lengthscale) {

  .gaussian_kernel_from_distances(X1, X2, lengthscale)
}


#------------------------------------------------------------------------------------------------------------------
#' Gaussian covariance kernel from optional precomputed squared distances
#'
#' @param X1 A numeric matrix.
#' @param X2 A numeric matrix.
#' @param lengthscale A numeric vector of Gaussian range parameters.
#' @param dist_sq Optional per-covariate squared distances from
#' `.gaussian_squared_distances()`.
#'
#' @return A covariance matrix.
#'
#' @keywords internal
.gaussian_kernel_from_distances <- function(X1, X2, lengthscale, dist_sq=NULL) {

  if(!is.null(dist_sq)){
    if(length(lengthscale) == 1){
      dist_scaled <- dist_sq / lengthscale
    }else{
      dist_scaled <- matrix(0, nrow = dim(dist_sq)[1], ncol = dim(dist_sq)[2])
      for(k in seq_along(lengthscale)){
        dist_scaled <- dist_scaled + dist_sq[, , k] / lengthscale[k]
      }
    }

    return(exp(-pmax(dist_scaled, 0)))
  }

  X1_scaled <- sweep(as.matrix(X1), 2, sqrt(lengthscale), "/")
  X2_scaled <- sweep(as.matrix(X2), 2, sqrt(lengthscale), "/")
  dist_sq <- outer(rowSums(X1_scaled^2), rowSums(X2_scaled^2), "+") -
    2 * tcrossprod(X1_scaled, X2_scaled)
  exp(-pmax(dist_sq, 0))
}


#------------------------------------------------------------------------------------------------------------------
#' Precompute per-covariate squared distances
#'
#' @param X1 A numeric matrix.
#' @param X2 A numeric matrix.
#'
#' @return A matrix for one covariate or an array with dimensions
#' nrow(X1) x nrow(X2) x ncol(X1) for multiple covariates.
#'
#' @keywords internal
.gaussian_squared_distances <- function(X1, X2) {

  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)

  if(ncol(X1) != ncol(X2)){
    stop("X1 and X2 must have the same number of columns.")
  }

  if(ncol(X1) == 1){
    out <- outer(X1[, 1], X2[, 1], "-")^2
    storage.mode(out) <- "double"
    return(out)
  }

  out <- array(0, dim = c(nrow(X1), nrow(X2), ncol(X1)))
  for(k in seq_len(ncol(X1))){
    out[, , k] <- outer(X1[, k], X2[, k], "-")^2
  }
  storage.mode(out) <- "double"
  out
}


#------------------------------------------------------------------------------------------------------------------
#' Cholesky factorization with adaptive jitter
#'
#' @param K A covariance matrix.
#' @param initial_jitter Initial diagonal jitter.
#' @param max_tries Maximum number of jitter increases.
#'
#' @return A list with the Cholesky factor and final jitter.
#'
#' @keywords internal
.chol_with_jitter <- function(K, initial_jitter=sqrt(.Machine$double.eps), max_tries=8) {

  jitter <- 0
  for(i in seq_len(max_tries)){
    chol_result <- tryCatch(chol(K + diag(jitter, nrow(K))),
                            error = function(cond) NULL)
    if(!is.null(chol_result)){
      return(list(chol = chol_result, jitter = jitter))
    }
    jitter <- if(i == 1) initial_jitter else jitter * 10
  }

  stop("Exact GP covariance matrix is not positive definite, even after adding jitter.")
}


#------------------------------------------------------------------------------------------------------------------
#' Standardize GP backend names
#'
#' @param gp_backend A character value indicating which GP backend to use.
#'
#' @return A normalized backend name.
#'
#' @keywords internal
.standardize_gp_backend <- function(gp_backend) {

  gp_backend <- match.arg(gp_backend, c("lagp", "exact_gaussian", "exact"))

  if(identical(gp_backend, "exact")){
    return("exact_gaussian")
  }

  gp_backend
}


#------------------------------------------------------------------------------------------------------------------
#' Fit a Gaussian process model to data
#'
#' This function fits a Gaussian process model to the provided data.
#'
#' @param X_data A data.frame containing the covariates used for model fitting. Should have at least timepoint, timepoint_n, and od columns.
#' @param delete A logical value. If TRUE, the Gaussian process model will be deleted after fitting.
#' @param gp_backend A character value indicating which GP backend to use.
#' @param gp_lengthscale Optional Gaussian-kernel lengthscale value or vector.
#' @param gp_nugget Optional observation-noise variance used as a nugget for
#' numerical stability.
#' @param gp_replicate_strategy A character value indicating whether replicate
#' rows with identical covariates should be used as raw observations or
#' aggregated to unbiased mean estimates with heteroskedastic nugget weights.
#' @param optimize_hyperparameters A logical value. If TRUE, estimate the
#' lengthscale and nugget by maximum likelihood.
#'
#' @return The input object, with the Gaussian process model information as a list
#'
#' @keywords internal
.fit_gp <- function(X_data, delete=FALSE, gp_backend="lagp", gp_lengthscale=NULL,
                    gp_nugget=NULL, gp_replicate_strategy="raw",
                    optimize_hyperparameters=TRUE) {

  gp_backend <- .standardize_gp_backend(gp_backend)

  if(identical(gp_replicate_strategy, "aggregate")){
    X_data <- .aggregate_replicate_timepoints(X_data)

    if(is.null(gp_nugget)){
      gp_nugget <- attr(X_data, "gp_nugget_weights")
    }
  }

  # Prepare X matrix
  X_matrix <- X_data %>%

    # In principle everything else should be timepoint and covariates
    select(-c(timepoint_n, od)) %>%
    as.matrix()

  # Gather od data for modelling
  y_vector <- X_data$od

  # Fit GP model
  fitted_gp <- switch(gp_backend,
                      lagp = .fit_lagp(X_matrix, y_vector, delete=FALSE,
                                       optimize_hyperparameters=optimize_hyperparameters),
                      exact_gaussian = .fit_exact_gaussian_gp(
                        X_matrix, y_vector,
                        lengthscale=gp_lengthscale,
                        nugget=gp_nugget,
                        optimize_hyperparameters=optimize_hyperparameters
                      ))

  # Return fitted model
  return(fitted_gp)
}


#------------------------------------------------------------------------------------------------------------------
#' Aggregate replicate measurements with identical GP covariates
#'
#' @param X_data A data.frame containing timepoint, timepoint_n, od, and
#' optional model covariates such as contrast_val.
#'
#' @return Aggregated data with nugget weights stored in the
#' gp_nugget_weights attribute.
#'
#' @keywords internal
.aggregate_replicate_timepoints <- function(X_data) {

  grouping_cols <- setdiff(colnames(X_data), c("curve_id", "timepoint_n", "od"))

  out <- X_data %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(timepoint_n = min(timepoint_n),
              od_mean = mean(od, na.rm = TRUE),
              n_replicates = sum(!is.na(od)),
              od_var = stats::var(od, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(timepoint_n)

  se2 <- out$od_var / pmax(out$n_replicates, 1)
  positive_se2 <- se2[is.finite(se2) & se2 > 0]

  if(length(positive_se2) == 0){
    nugget_weights <- rep(1, nrow(out))
  }else{
    fallback_se2 <- stats::median(positive_se2)
    se2[!is.finite(se2) | se2 <= 0] <- fallback_se2
    nugget_weights <- se2 / fallback_se2
  }

  out <- out %>%
    transmute(across(all_of(grouping_cols)),
              timepoint_n = timepoint_n,
              od = .data$od_mean)

  attr(out, "gp_nugget_weights") <- pmax(nugget_weights, sqrt(.Machine$double.eps))
  out
}


#------------------------------------------------------------------------------------------------------------------
#' Generate predictions from a Gaussian Process model
#'
#' This function generates predictions from a Gaussian Process model fitted to an object.
#' The function supports "alternative" and "null" models. The "alternative" model predicts
#' based on both time and contrast, while the "null" model predicts based on time only.
#'
#' @param XX The covariate matrix for which we want ot make predictions
#' @param lagp_fit The output of an laGP fit.
#' @param complete_sigma A logical value indicating whether the full covariance matrix should be returned. This should be set to TRUE if sampling the posterior is desired.
#' @param delete A logical value. If TRUE, the Gaussian Process model will be deleted from the object after prediction.
#'
#'
#' @return The mean predicted values.
#'
#' @seealso [laGP::predGPsep()]
#'
#' @keywords internal
.predict_lagp <- function(XX, lagp_fit, complete_sigma=FALSE, delete=TRUE){

  # At this point we can directly fit
  predicted_lagp <- laGP::predGPsep(gpsepi = lagp_fit$gpsepi, XX = XX, lite = !complete_sigma, nonug = complete_sigma)

  # Update laGP object
  lagp_fit$prediction <- predicted_lagp

  if (delete == TRUE) {
    laGP::deleteGPsep(lagp_fit$gpsepi)
    lagp_fit$gpsepi <- NULL
  }

  return(lagp_fit)
}


#------------------------------------------------------------------------------------------------------------------
#' Generate predictions from a Gaussian Process model
#'
#' This function generates predictions from a Gaussian Process model fitted to an object.
#' The function supports "alternative" and "null" models. The "alternative" model predicts
#' based on both time and contrast, while the "null" model predicts based on time only.
#'
#' @param XX The covariate matrix for which we want ot make predictions
#' @param lagp_fit The output of an laGP fit.
#' @param complete_sigma A logical value indicating whether the full covariance matrix should be returned. This should be set to TRUE if sampling the posterior is desired.
#' @param estimate_derivatives A logical value indicating if the empirical first and second derivatives should be returned.
#' @param delete A logical value. If TRUE, the Gaussian Process model will be deleted from the object after prediction.
#'
#'
#' @return The mean predicted values.
#'
#' @seealso [laGP::predGPsep()]
#'
#' @keywords internal
.predict_gp <- function(XX, lagp_fit, complete_sigma=FALSE, prepare_dataframe=FALSE, estimate_derivatives = TRUE, delete=TRUE){

   # At this point we can directly fit
   gp_backend <- if(is.null(lagp_fit$backend)) "lagp" else lagp_fit$backend
   predicted_lagp <- switch(gp_backend,
                            exact_gaussian = .predict_exact_gaussian_gp(XX, lagp_fit,
                                                                        complete_sigma=complete_sigma),
                            lagp = .predict_lagp(XX, lagp_fit,
                                                 complete_sigma=complete_sigma,
                                                 delete=delete))

   # Update laGP object
   if(prepare_dataframe){

     predicted_lagp$prediction$prediction_dataframe <- XX %>%
       as.data.frame() %>%

       mutate(mean = predicted_lagp$prediction$mean, # Mean of posterior

              diag.sigma = ifelse(complete_sigma,
                                  diag(predicted_lagp$prediction$Sigma),
                                  predicted_lagp$prediction$s2), # Gather diagonal of covariance
       )

     # check if derivatives should be estimated
     if(estimate_derivatives){
       if(identical(gp_backend, "exact_gaussian")){
         analytic_derivatives <- .exact_gaussian_derivatives(XX, predicted_lagp)
         predicted_lagp$prediction$prediction_dataframe <- bind_cols(
           predicted_lagp$prediction$prediction_dataframe,
           analytic_derivatives
         )
       }else{
         predicted_lagp$prediction$prediction_dataframe <- .estimate_derivatives(predicted_lagp$prediction$prediction_dataframe)
       }
     }

   }


   return(predicted_lagp)

}


#------------------------------------------------------------------------------------------------------------------
#' Sample the posterior GP model for new curves
#'
#' @param lagp_fit The output of an laGP fit.
#' @param n_samples A numberic value indicating the number of curves to sample from posterior
#'
#'
#' @return A data.frame with sample_id, timepoint, and mean columns.
#'
#'
#' @seealso [laGP::predGPsep()]
#' @importFrom mvtnorm rmvt
#' @keywords internal
.sample_posterior_fit <- function(lagp_fit, n_samples=100){

   # Sample according to t-distribution
   ZZ <- rmvt(n_samples, lagp_fit$prediction$Sigma, lagp_fit$prediction$df)
   ZZ <- ZZ + t(matrix(rep(lagp_fit$prediction$mean, n_samples), ncol = n_samples))

   # Makes things a bit easier
   ZZ <- t(ZZ)

   # Label samples
  colnames(ZZ) <- paste0("posterior_sample_", 1:n_samples)

   # Now format nicely
   samples_df <- ZZ %>%
     as.data.frame() %>%
     mutate(timepoint = lagp_fit$prediction$prediction_dataframe[, "timepoint"]) %>%

     pivot_longer(cols=-timepoint, names_to = "sampled_id", values_to = "mean")
   return(samples_df)


}



#------------------------------------------------------------------------------------------------------------------
#' Fit at GP model and predict the mean process
#'
#' This function fits a Gaussian Process model to a data.frame and then generates predictions from the fitted model.
#'
#' @param df A data.frame with the data that is to be modelled. Must have timepoint, timepoint_n, curve_id, and od fields.
#' @param t_steps the number of timepoints within the range od df$timepoints to predict
#' @param complete_sigma A logical indicating whether the complete covariance matrix should be returned. This is necessary if you want to re-sample the posterior.
#' @param prepare_dataframe A logical indicating whether a data.frame with the mean data should be prepared.
#' @param estimate_derivatives A logical indicating whether the empirical derivatives should be computed.
#' @param delete A logical indicating whether the GP model should be deleted from memory FALSE.
#' @param pred_for_t A numeric value inidicating a specific timepoint t for which OD should be predicted
#' @param gp_backend A character value indicating which GP backend to use. The
#' default "lagp" preserves previous DGrowthR behavior. The "exact_gaussian"
#' backend uses dense Cholesky GP inference with a Gaussian kernel, marginal
#' log-likelihood, explicit nugget, and analytic derivatives for short
#' trajectories. The old "exact" value is accepted as an alias.
#' @param gp_lengthscale Optional Gaussian-kernel lengthscale value or vector.
#' @param gp_nugget Optional observation-noise variance used as a nugget for
#' numerical stability.
#' @param gp_replicate_strategy A character value indicating whether replicate
#' rows with identical covariates should be used as raw observations or
#' aggregated to unbiased mean estimates with heteroskedastic nugget weights.
#' @param optimize_hyperparameters A logical value. If TRUE, estimate the
#' lengthscale and nugget by maximum likelihood. Set to FALSE for faster
#' approximate fits that use laGP's initial hyperparameter values.
#'
#' @return An object with the fitted model and prediction results. The results include the predictions (p) and the prediction inputs (XX).
#'
#' @rdname fit_predict_gp
#' @export
fit_predict_gp <- function(df, t_steps=1000, complete_sigma=FALSE, prepare_dataframe=FALSE, estimate_derivatives=FALSE, delete=TRUE, pred_for_t=NULL, gp_backend=c("lagp", "exact_gaussian", "exact"), gp_lengthscale=NULL, gp_nugget=NULL, gp_replicate_strategy=c("raw", "aggregate"), optimize_hyperparameters=TRUE){

  gp_backend <- .standardize_gp_backend(gp_backend)
  gp_replicate_strategy <- match.arg(gp_replicate_strategy)

  # Check if curve_id field is present
  if("curve_id" %in% colnames(df)){
    df <- df %>%
      select(-curve_id)
  }

  if("contrast_val" %in% colnames(df)){
    contrast_value <- unique(df$contrast_val)
  }

  # Fit the object
  gpobj <- .fit_gp(df, delete = FALSE,
                   gp_backend=gp_backend,
                   gp_lengthscale=gp_lengthscale,
                   gp_nugget=gp_nugget,
                   gp_replicate_strategy=gp_replicate_strategy,
                   optimize_hyperparameters=optimize_hyperparameters)


  # Determine if a simple or complex matrix was supplied
  if("contrast_val" %in% colnames(df)){
    time_vector_double <- rep(seq(from=min(gpobj$X[, "timepoint"]),
                           to=max(gpobj$X[, "timepoint"]),
                           length.out=t_steps), 2)

    # Prepare an XX for prediction
    XX <- data.frame(timepoint = time_vector_double,
                     contrast_val = rep(c(0,1), each=t_steps)) %>%
      as.matrix()

  }else{

    # Prepare an XX for prediction
    XX <- data.frame(timepoint = seq(from=min(gpobj$X[, "timepoint"]),
                                     to=max(gpobj$X[, "timepoint"]),
                                     length.out=t_steps))


    # Include specific timepoint in the prediction if requested
    if(!is.null(pred_for_t)){

      new_timepoints <- sort(unique(c(pred_for_t, XX$timepoint)))

      XX <- data.frame(timepoint = new_timepoints)

    }
  }


  # Predict mean function
  gpobj <- .predict_gp(XX=XX,
                       lagp_fit=gpobj,
                       complete_sigma = complete_sigma,
                       prepare_dataframe=prepare_dataframe,
                       estimate_derivatives=estimate_derivatives,
                       delete = delete)

  return(gpobj)
}
