\name{fit_gp}
\alias{fit_gp}
\alias{predict_gp}
\alias{assess_gp_fit}
\title{Fit, predict, and assess Gaussian processes}
\usage{
fit_gp(object, model_covariate = "curve_id",
  downsample_every_n_timepoints = 1,
  gp_backend = c("lagp", "exact_gaussian", "exact"),
  gp_lengthscale = NULL, gp_nugget = NULL,
  gp_replicate_strategy = c("raw", "aggregate"),
  optimize_hyperparameters = TRUE, n_cores = 1)

predict_gp(object, predict_n_steps = 100, prediction_timepoints = NULL,
  complete_sigma = FALSE, estimate_derivatives = TRUE)

assess_gp_fit(object, level = 0.95,
  interval_type = c("predictive", "latent"),
  review_duration = 0.5, review_points = 4L,
  fail_duration = 1, fail_excess = 1)
}
\arguments{
\item{object}{A DGrowthR object.}
\item{model_covariate}{Metadata column defining GP groups.}
\item{downsample_every_n_timepoints}{Use every nth observed timepoint.}
\item{gp_backend}{The laGP or exact Gaussian backend.}
\item{gp_lengthscale}{Optional exact-GP lengthscale.}
\item{gp_nugget}{Optional exact-GP nugget.}
\item{gp_replicate_strategy}{Use raw or aggregated replicates.}
\item{optimize_hyperparameters}{Optimize lengthscale and nugget.}
\item{n_cores}{Number of fitting workers.}
\item{predict_n_steps}{Number of equally spaced prediction timepoints.}
\item{prediction_timepoints}{Optional explicit numeric prediction timepoints.}
\item{complete_sigma}{Store the full latent posterior covariance.}
\item{estimate_derivatives}{Calculate first and second derivatives.}
\item{level}{Two-sided interval probability.}
\item{interval_type}{Predictive interval for group-mean observations or latent interval for the mean process.}
\item{review_duration}{Outside-run duration threshold for review, in timepoint units.}
\item{review_points}{Outside-run point-count threshold for review.}
\item{fail_duration}{Outside-run duration threshold for failure, in timepoint units.}
\item{fail_excess}{Maximum relative interval-excess threshold for failure.}
}
\description{
\code{fit_gp()} stores serializable group fits without live laGP handles.
\code{predict_gp()} reuses those fits and stores prediction data.
\code{assess_gp_fit()} evaluates group-mean observations at their original
timepoints and reports residual, coverage, AUC, and sustained outside-band
diagnostics.
}
\value{
\code{fit_gp()} and \code{predict_gp()} return a DGrowthR object;
\code{assess_gp_fit()} returns a list of summary, runs, and points data frames.
}
