#' @title Estimate growth parameters
#'
#' @description This function estimates a collection of growth parameters from the mean GP fits. If necessary, it also runs the GP fits.
#'
#' @param object A DGrowthR object containing preprocessed data.
#' @param model_covariate A string indicating a covariate in metadata to pool growth curves for GP modelling.
#' @param predict_n_steps A numeric value indicating the number of timepoints to sample from the mean posterior of the GP fit.
#' @param downsample_every_n_timepoints A numeric value indicating that the OD from every n timepoint should be used for GP fit. Might seep up fitting.
#' @param sample_posterior_gpfit A logical value. Indicates if the posterior GP should be sampled for new growth curves and parameters are estimated
#' from these sampled curves. Allows to estimate mean and standard deviation of growth parameters. Really only makes sense if model_covariate pools more
#' than one growth curve per group.
#' @param od_auc_at_t A numeric value indicating a specific timepoint for which the predicted OD and AUC should be returned.
#' @param sample_n_curves A numeric values. If sample_posterior_gpfit is TRUE, then sample_n_curves are sampled from posterior.
#' @param save_gp_data A logical value, indicating if the mean GP values and GP fit parameters should be saved to object.
#' @param n_cores A numeric values indicating the number of cores the user wants to use to model curves in parallel.
#' @param detect_polyauxic A logical value indicating if the function should attempt to detect polyauxic growth curves. This is done by looking for multiple peaks in the second derivative of the mean GP fit. If TRUE, then additional columns to the growth parameter table will be added with the relevant information.
#' @param polyauxic_ratio_threshold A numeric value. Potential growth phases are assessed by comparing their respective carrying capacity. If the carrying capacity of a potential second phase is at least this fraction of the carrying capacity of the first phase, then the second phase is considered a true growth phase. Only relevant if detect_polyauxic is TRUE.
#' @param polyauxic_parameter_criteria A character value. Either "growth_rate" or "carrying_capacity". Indicates which parameter should be used to assess the relevance of potential growth phases. Only relevant if detect_polyauxic is TRUE.
#' @return Updated DGrowthR object with updated growth_parameters and, if requested, gpfit_info slots.
#'
#'
#' @export
setGeneric("estimate_growth_parameters", function(object,
                                                  model_covariate="curve_id",
                                                  predict_n_steps=100,
                                                  downsample_every_n_timepoints=1,
                                                  sample_posterior_gpfit=FALSE,
                                                  od_auc_at_t = NULL,
                                                  sample_n_curves=100,
                                                  save_gp_data = FALSE,
                                                  n_cores=1,
                                                  detect_polyauxic = FALSE,
                                                  polyauxic_ratio_threshold = 0.2,
                                                  polyauxic_parameter_criteria = "carrying_capacity"){
  standardGeneric("estimate_growth_parameters")
})

#' @rdname estimate_growth_parameters
#' @export
setMethod(
  f = "estimate_growth_parameters",
  signature = "DGrowthR",
  definition = function(object, model_covariate="curve_id", predict_n_steps=100,
                        downsample_every_n_timepoints=1,
                        sample_posterior_gpfit=FALSE,
                        od_auc_at_t = NULL,
                        sample_n_curves=100, save_gp_data = FALSE, n_cores=1,
                        detect_polyauxic = FALSE, polyauxic_ratio_threshold = 0.2,
                        polyauxic_parameter_criteria = "carrying_capacity"){


    # Some checks first
    if(!object@preprocessed){
      warning("Growth data has not been pre-processed by DGrowthR. This may lead to some problems. Try running preprocess_data() first")
    }

    if(!object@log_od){
      warning("OD data is not log-transformed. This may lead to innacurate estimates. Try running preprocess_data(..., log_transform=TRUE)")
    }

    # Determine if we need to run (or re-run) the GP fit due to not fit data found or re-sampling requested.
    if(length(object@gpfit_info) == 0 || sample_posterior_gpfit){
      message("[DGrowthR::estimate_growth_parameters] >> Fitting GP models and estimating growth parameters...")
      object <- .fit_and_estimate(object,  model_covariate, predict_n_steps, downsample_every_n_timepoints, sample_posterior_gpfit, od_auc_at_t, sample_n_curves, save_gp_data, n_cores, detect_polyauxic, polyauxic_ratio_threshold, polyauxic_parameter_criteria)

    # Check if paramters for GP fit have changed
      }else if(object@gpfit_info$gpfit_parameters$model_covariate != model_covariate || object@gpfit_info$gpfit_parameters$predict_n_steps != predict_n_steps || object@gpfit_info$gpfit_parameters$downsample_every_n_timepoints != downsample_every_n_timepoints || object@gpfit_info$gpfit_parameters$detect_polyauxic != detect_polyauxic || object@gpfit_info$gpfit_parameters$polyauxic_ratio_threshold != polyauxic_ratio_threshold || object@gpfit_info$gpfit_parameters$polyauxic_parameter_criteria != polyauxic_parameter_criteria){

      message("[DGrowthR::estimate_growth_parameters] >> Fitting new GP models due to change in hyperparameters...")
      object <- .fit_and_estimate(object,  model_covariate, predict_n_steps, downsample_every_n_timepoints, sample_posterior_gpfit, od_auc_at_t, sample_n_curves, save_gp_data, n_cores, detect_polyauxic, polyauxic_ratio_threshold, polyauxic_parameter_criteria)

    # If there is GP data, no changes to the modelling are requested, and no sampling posterior requested, then use available gpfits
    }else{
      message("[DGrowthR::estimate_growth_parameters] >> Using stored GP data...")

      growth_parameters <- object@gpfit_info$gpfit_data %>%
        group_by(gpfit_id) %>%
        group_modify(~.characterize_growth(.x)) %>%
        ungroup()


      slot(object, "growth_parameters") <- growth_parameters

      }

  message("[DGrowthR::estimate_growth_parameters] >> Finished!")
  return(object)

  }
)

#------------------------------------------------------------------------------------------------------------------
#' Fit GP and gather growth parameters
#'
#' @param object A DGrowthR object containing preprocessed data.
#' @param model_covariate A string indicating a covariate in metadata to pool growth curves for GP modelling.
#' @param predict_n_steps A numeric value indicating the number of timepoints to sample from the mean posterior of the GP fit.
#' @param downsample_every_n_timepoints A numeric value indicating that the OD from every n timepoint should be used for GP fit. Might seep up fitting.
#' @param sample_posterior_gpfit A logical value. Indicates if the posterior GP should be sampled for new growth curves and parameters are estimated
#' from these sampled curves. Allows to estimate mean and standard deviation of growth parameters. Really only makes sense if model_covariate pools more
#' than one growth curve per group.
#' @param od_auc_at_t A numeric value indicating a specific timepoint for which the predicted OD and AUC should be returned.
#' @param sample_n_curves A numeric values. If sample_posterior_gpfit is TRUE, then sample_n_curves are sampled from posterior.
#' @param save_gp_data A logical value, indicating if the mean GP values and GP fit parameters should be saved to object.
#' @param n_cores A numeric values indicating the number of cores the user wants to use to model curves in parallel.
#' @param detect_polyauxic A logical value indicating if the function should attempt to detect polyauxic growth curves. This is done by looking for multiple peaks in the second derivative of the mean GP fit. If TRUE, then and additional dataframe with the timepoints of the beginning and end of each phase is included. If set to TRUE, then save_gp_data must also be set to TRUE.
#' @param polyauxic_ratio_threshold A numeric value. Potential growth phases are assessed by comparing their respective carrying capacity. If the carrying capacity of a potential second phase is at least this fraction of the carrying capacity of the first phase, then the second phase is considered a true growth phase. Only relevant if detect_polyauxic is TRUE.
#' @param polyauxic_parameter_criteria A character value. Either "growth_rate" or "carrying_capacity". Indicates which parameter should be used to assess the relevance of potential growth phases. Only relevant if detect_polyauxic is TRUE.
#'
#' @return Updated DGrowthR object
#'
#' @seealso [laGP::predGPsep()]
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom magrittr %>%
#'
#' @keywords internal
.fit_and_estimate <- function(object,  model_covariate, predict_n_steps,
                              downsample_every_n_timepoints,
                              sample_posterior_gpfit, od_auc_at_t,
                              sample_n_curves, save_gp_data, n_cores,
                              detect_polyauxic, polyauxic_ratio_threshold, polyauxic_parameter_criteria){

  # Gather od and metadata
  #complete_od <- object@od_data
  metadata <- object@metadata

  # Check that covariate is in metadata
  if(!model_covariate %in% colnames(metadata)){

    stop(paste("Covariate", model_covariate, "not found in metadata."))

  }

  # Determine the groups we are trying to model
  if(model_covariate == "curve_id"){

    gpfit_metadata <- metadata %>%
      mutate(gpfit_id = curve_id) %>%
      select(curve_id, gpfit_id)

  }else{

    gpfit_metadata <- metadata %>%
      rename("gpfit_id" = all_of(model_covariate)) %>%
      select(curve_id, gpfit_id)
  }

  # Gather a vector of all of the groups we want to model
  model_groups <- metadata %>%
    select(all_of(model_covariate)) %>%
    distinct() %>%
    unlist()

  # Print how many models we will create
  message(paste("[DGrowthR::estimate_growth_parameters] >> Modelling the", model_covariate, "field from metadata.", length(model_groups), "models will be created."))

  # Registering parellelization
  #doParallel::registerDoParallel(n_cores)

  cl <- snow::makeSOCKcluster(n_cores)
  doSNOW::registerDoSNOW(cl)

  # Set up progress bar
  pb <- txtProgressBar(max = length(model_groups), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)

  # Iterate over groupings and fit GP and estimate growth params
  estimate_list <- foreach::foreach(i = model_groups, .combine = c,
                                    .multicombine = TRUE, .init = c(),
                                    .options.snow = list(progress = progress)) %dopar% {

    # Gather the relevant growth data
    od_for_gp <- gather_od_data(object, model_covariate, i, downsample_every_n_timepoints) %>%
      select(timepoint, timepoint_n, od)


    # Check if curve is flat
    if(.check_flat(od_for_gp)){
      message(paste("[DGrowthR::estimate_growth_parameters] >> flat curve detected for", i, "..."))
      status <- "flat_curve"

      growth_params <- data.frame(gpfit_id = i,
                                  gpfit_converged = status)
      # Prepare output
      out_list <- list("growth_parameters" = growth_params)
      return(out_list)
    }

    # Fit the GP
    gpfit <- tryCatch(
      {

        fit_predict_gp(od_for_gp,
                       t_steps=predict_n_steps,
                       complete_sigma=sample_posterior_gpfit,
                       prepare_dataframe=TRUE,
                       estimate_derivatives=TRUE,
                       delete=TRUE,
                       pred_for_t=od_auc_at_t)

        },
      error = function(cond){
        print(i)
      })

    # Determine if growth parameters should be determined by resampling posterior or by mean curve
    if(sample_posterior_gpfit){

      # Sample from posterior fit
      sampled_mean_gps <- .sample_posterior_fit(gpfit, n_samples=sample_n_curves)

      # Gather the derivatives for each sampled curve
      sampled_mean_gps <- sampled_mean_gps %>%
        group_by(sampled_id) %>%
        group_modify(~.estimate_derivatives(.x)) %>%
        ungroup()

      # Characterize growth for each sampled curve
      growth_params <- sampled_mean_gps %>%
        group_by(sampled_id) %>%
        group_modify(~.characterize_growth(.x)) %>%

        # Gather summary statistics
        ungroup() %>%
        select(-sampled_id) %>%
        summarise(across(.cols=where(is.numeric),
                  list(mean = ~mean(.x, na.rm=TRUE), sd = ~sd(.x, na.rm = TRUE)), # gather average and standard dev,
                  .names = "{.col}.{.fn}"))

    }else{
      # Else, the assumption is that we provide the mean OD
      growth_params <- .characterize_growth(gpfit$prediction$prediction_dataframe,
                                            od_auc_at_t,
                                            detect_polyauxic,
                                            polyauxic_ratio_threshold,
                                            polyauxic_parameter_criteria)
    }


    # Add identifying column
    growth_params <- growth_params %>%
      mutate(gpfit_id = i,
             gpfit_converged = ifelse(gpfit$mle$conv == 0, "yes", "no")) %>%
      relocate(gpfit_id)


    # Prepare output
    out_list <- list("growth_parameters" = growth_params)

    if(save_gp_data){
      out_list$gp_data <- gpfit$prediction$prediction_dataframe
      out_list$gp_data$gpfit_id <- i
    }

    return(out_list)
  }

  close(pb)
  snow::stopCluster(cl)

  # Update object with growth parameters
  message("[DGrowthR::estimate_growth_parameters] >> Updating growth_parameters slot...")
  #print(length(estimate_list))
  if(length(model_groups) == 1){

    slot(object, "growth_parameters") <- estimate_list$growth_parameters

  }else{

    slot(object, "growth_parameters") <- bind_rows(estimate_list[names(estimate_list) == "growth_parameters"]) %>%
      relocate(gpfit_converged, .after = last_col())

  }


  # Update with gp if necessary
  if(save_gp_data){
    message("[DGrowthR::estimate_growth_parameters] >> Updating gpfit_info slot...")

    gpfit_list <- list("gpfit_parameters" = list("model_covariate"=model_covariate,
                                                "predict_n_steps"=predict_n_steps,
                                                "downsample_every_n_timepoints"=downsample_every_n_timepoints))

    # Update object with growth parameters
    if(length(model_groups) == 1){

      gpfit_list$gpfit_data <- estimate_list$gp_data

    }else{

      gpfit_list$gpfit_data <- bind_rows(estimate_list[names(estimate_list) == "gp_data"])

    }

    slot(object, "gpfit_info") <- gpfit_list
  }


  return(object)

}


#------------------------------------------------------------------------------------------------------------------
#' Estimate growth parameters
#'
#' @param df A data.frame with timepoint and od fields.
#'
#' @return A logical value indicating if the curve is flat
#'
#'
#' @keywords internal
.check_flat <- function(df){

  nunique_od_vals <- df %>%
    select(timepoint, od) %>%
    distinct() %>%
    select(od) %>%
    n_distinct()

  # Make sure there are more than 3 unique values
  return(nunique_od_vals <= 3)

}


#------------------------------------------------------------------------------------------------------------------
#' Estimate growth parameters
#'
#' @param df A data.frame with timepoint and mean fields.
#' @param od_auc_at_t A numeric value for which od and AUC are calculated for
#' @param detect_polyauxic A logical value indicating if the function should attempt to detect polyauxic growth curves. This is done by looking for multiple peaks in the second derivative of the mean GP fit. If TRUE, then and additional dataframe with the timepoints of the beginning and end of each phase is included.
#' @param polyauxic_ratio_threshold A numeric value. Potential growth phases are assessed by comparing their respective carrying capacity. If the carrying capacity of a potential second phase is at least this fraction
#' @param polyauxic_parameter_criteria A character value. Either "growth_rate" or "carrying_capacity". Indicates which parameter should be used to assess the relevance of potential growth phases. Only relevant if detect_polyauxic is TRUE.
#'
#' @return A data.frame with the requested growth parameters for the curve given as input
#'
#'
#' @keywords internal
.characterize_growth <- function(df, od_auc_at_t=NULL, detect_polyauxic = FALSE,
                                 polyauxic_ratio_threshold = 0.2,
                                 polyauxic_parameter_criteria = "carrying_capacity"){

  # Maximum od and max od timepoint
  max_od_data <- .determine_max_growth(df)

  # AUC
  auc_od_data <- .determine_auc(df)

  # Growth loss
  growth_loss <- .determine_growth_loss(df, max_od_data)

  # Determine growth rate and death rate
  max_growth_rate <- .determine_max_growth_rate(df, max_od_data)

  # Determine growth rate and death rate
  max_death_rate <- .determine_max_death_rate(df, max_od_data)

  # Determining the end of lag-phase and the possible start of stationary phase
  # As a heuristic, we look for the the moment of fastest aceleration before the maximum growth rate
  # and the moment of greatest deceleration after the maximum growth rate.
  end_lagphase <- .determine_lag_phase(df, max_growth_rate)
  start_statphase <- .determine_stationary_phase(df, max_growth_rate, max_od_data)


  # Join all growth parameters
  gparams <- bind_cols(list(max_od_data, auc_od_data, growth_loss,
                         max_growth_rate, max_death_rate, end_lagphase,
                         start_statphase))

  # If data is requested for a specific timepoint
  if(!is.null(od_auc_at_t)){

    # Rename od
    od_colname <- paste0("OD_", od_auc_at_t)

    # Auc until timepoint
    auc_until_t <- .determine_auc(df, od_auc_at_t)

    # OD at timepoint
    od_at_t <- df %>%
      filter(timepoint == od_auc_at_t) %>%
      select(mean)
    colnames(od_at_t) <- od_colname


    # Add to gparams
    gparams <- bind_cols(list(gparams, auc_until_t, od_at_t))

  }

  # If requested, attempt to detect polyauxic growth curves
  if(detect_polyauxic){
    # Check if there are multiple peaks in the second derivative
    # If there are, then we have potential multiple growth phases. We can use the same functions as before to determine the growth parameters of the second phase, but we need to make sure that the carrying capacity of the second phase is at least a certain fraction of the first phase (as determined by polyauxic_ratio_threshold) to avoid overfitting.

    polyauxic_df <- .detect_polyauxic(df, polyauxic_ratio_threshold, polyauxic_parameter_criteria)
    gparams <- bind_cols(list(gparams, polyauxic_df))

  }

  return(gparams)
}

##------------------------------------------------------------------------------------------------------------------
#' Detect polyauxic growth curves
#' This function is heavily based on the implementation available on AMiGA.
#'
#' @param df A data.frame with timepoint, mean, first and second derivative fields.
#' @param thresh A numeric value indicating the minimum fraction of the carrying capacity of the first phase that the second phase must have to be considered a true growth phase.
#' @param varb_detect A character value. Either "growth_rate" or "carrying_capacity". Indicates which parameter should be used to assess the relevance of potential growth phases. Only relevant if detect_polyauxic is TRUE.
#'
#' @return A data frame.
#'
#' @keywords internal
.detect_polyauxic <- function(df, thresh, varb_detect){

  # Gather the mean, first and second derivatives as vectors for easier handling
  y <- df$mean
  d1 <- df %>% filter(!is.na(second_derivative)) %>% pull(first_derivative)
  d2 <- df %>% filter(!is.na(second_derivative)) %>% pull(second_derivative)
  time_vec <- df$timepoint

  # Step 1: Identify moments where the roots of the second derivative change
  # from positive to negative. These are potential inflection points that could
  # indicate the end of a growth phase and the start of a new one.
  potential_inflection_points <- which(diff(sign(d2)) != 0)


  # Step 2.a: Determine the inflection direction at each potential inflection point.
  inflection_directions <- sapply(potential_inflection_points, function(ii) {
    # Check if we are at least 1 spot away from the absolute end
    if (ii < (length(d2) - 1)) {
      return(sign(d2[ii + 1])) # Look forward
    } else {
      return(-1 * sign(d2[ii - 1])) # Look backward (edge case)
    }
  })

  # Speical case. If there are no inflection points we can return early with no growth phases detected.
  if(length(potential_inflection_points) == 0){
    return(data.frame(n_growth_phases = 1))
  }

  # Step 2.b: Pad the inflection points and directions to ensure that we have a
  # complete set of inflection points and directions for the entire curve.
  # Pad the left edge (Start)
  res_left <- .pad_inflections(potential_inflection_points, inflection_directions,
                  edge = 1, len = length(d2))

  potential_inflection_points <- res_left$ips
  inflection_directions <- res_left$its

  # Pad the right edge (End)
  res_right <- .pad_inflections(potential_inflection_points, inflection_directions,
                   edge = -1, len = length(d2))
  potential_inflection_points <- res_right$ips
  inflection_directions <- res_right$its

  # Step 3: Define bounds for each potential growth phase based on the inflection
  # points and directions.

  # Find indices where inflection direction is positive (potential start of growth phase)
  starts <- which(inflection_directions == 1)
  starts <- head(starts, -1) # Remove the last one to avoid out-of-bounds

  # Define the stops corresponding to each stop
  # We assume a pattern of Start, Peak, Stop, and we are only looking for the start of the next phase,
  # so we can just look for the next positive inflection point after the current start.
  # This is a simplification that assumes that there are no nested growth phases, but it should work for most cases.
  stops <- starts + 2

  # Step 4: Initialize vectors to store the start and stop timepoints of each
  # growth phase and their corresponding growht parameters
  n <- length(starts)
  t_left_vec <- numeric(n)
  t_right_vec <- numeric(n)
  K_vec <- numeric(n)
  r_vec <- numeric(n)
  r_left_vec <- numeric(n) # Growth rate at the left inflection point (start of phase)
  r_right_vec <- numeric(n) # Growth rate at the right inflection point (end of phase)
  attraction_vec <- numeric(n)

  # Iterate through each growth stage
  for (i in 1:n) {

    # Step 4.a: Get the indices from the 'potential_inflection_points' vector
    # 'starts' and 'stops' are indices pointing into 'potential_inflection_points'
    idx_L <- potential_inflection_points[starts[i]]
    idx_R <- potential_inflection_points[stops[i]]

    t_left_vec[i]  <- idx_L
    t_right_vec[i] <- idx_R

    # Step 4.b: Define the slice range for calculation
    # We look at the data BETWEEN the inflection points.
    # Start at L+1 to exclude the inflection point itself (start of phase)
    # End at R (include the end inflection point)
    range_idxs <- (idx_L + 1):idx_R

    # Safety check: ensure range is valid
    if (length(range_idxs) > 0) {

      # Total change in OD (Max value - Start value)
      current_y_slice <- y[range_idxs]
      start_y_val     <- y[range_idxs[1]] # Value at the start of the slice
      K_vec[i]         <- max(current_y_slice - start_y_val)

      # Max growth rate in this window
      current_d1_slice <- d1[range_idxs]
      r_vec[i]         <- max(current_d1_slice)

      #Growth rates at the bounds
      r_left_vec[i]  <- d1[idx_L]
      r_right_vec[i] <- d1[idx_R]

      # Attraction (Skew)
      # "Is the drop to the right bigger than the drop to the left?"
      dist_right <- abs(r_right_vec[i] - r_vec[i])
      dist_left  <- abs(r_left_vec[i]  - r_vec[i])

      if (dist_right > dist_left) {
        attraction_vec[i] <- -1
      } else {
        attraction_vec[i] <- 1
      }
    }
  }

  # Create the final DataFrame
  ret <- data.frame(
    t_left     = t_left_vec,
    t_right    = t_right_vec,
    K          = K_vec,
    r          = r_vec,
    r_left     = r_left_vec,
    r_right    = r_right_vec,
    attraction = attraction_vec
  )

  # Determine the primary and secondary parameter for growth phase assessment based on user input
  pvarb <- ifelse(varb_detect == "carrying_capacity", "K", "r")
  svarb <- ifelse(varb_detect == "carrying_capacity", "r", "K")
  if(any(ret[[pvarb]] > 0)){
    # Step 5: Iteratively detect phases and merge if necessary based on the ratio threshold and the selected parameter criteria (growth rate or carrying capacity)

    ret <- .detect_phases(ret, ratio_threshold = thresh, varb=pvarb, second_varb=svarb)
  }else {
    # If no positive growth, merge everything into one
    while (nrow(ret) > 1) {
      ret <- .merge_phases(ret, 1, 2)
    }
  }

  # Convert Indices to Time
  # equivalent to: ret.iloc[:,0].apply(lambda i: x[int(i)])
  ret$t_left  <- time_vec[as.integer(ret$t_left)]
  ret$t_right <- time_vec[as.integer(ret$t_right)]

  # Drop attraction column
  ret$attraction <- NULL

  # Pivot wider
  ret_onerow <- ret %>%
    rename("t_start" = "t_left",
           "t_end" = "t_right",
           "max_growth" = "K",
           "max_growth_rate" = "r",
           "growth_rate_start" = "r_left",
           "growth_rate_end" = "r_right") %>%

    mutate(phase = paste0("phase", row_number())) %>%
    pivot_wider(names_from = phase,
                values_from = -phase,
                names_glue = "{.value}.{phase}") %>%

    mutate(n_growth_phases = nrow(ret)) %>%
    relocate(n_growth_phases)

  return(ret_onerow)

}


##------------------------------------------------------------------------------------------------------------------
#' Iteratively detect phases
#' This function is heavily based on the implementation available on AMiGA.
#'
#' @param df A data.frame with the parameters of each potential growth phase.
#' @param ratio_threshold A numeric value indicating the minimum fraction of the carrying capacity of the first phase that the second phase must have to be considered a true growth phase.
#' @return A data frame with the detected growth phases.
#' @keywords internal
.detect_phases <- function(df, ratio_threshold, varb="K", second_varb="r"){

  # While the smallest phase is smaller than threshold * max phase
  while (min(df[[varb]]) < ratio_threshold * max(df[[varb]])) {

    # 1. Force edge attraction
    # First row: Attract forward (1)
    # Last row: Attract backward (-1)
    df <- df[order(df$t_left), ]
    df$attraction[1] <- 1
    df$attraction[nrow(df)] <- -1

    # 2. Find the "weakest" phase
    # Sort by K (varb) then r (second_varb) to find the smallest one
    # We need to preserve original row numbers to know who is who,
    # so we add a temporary ID column.
    df$temp_id <- 1:nrow(df)

    sorted_ret <- df[order(df[[varb]], df[[second_varb]]), ]

    # Get the row index (position) of the smallest item
    idx_to_merge <- sorted_ret$temp_id[1]

    # 3. Determine target
    # Look at the attraction value of that specific row
    att <- df$attraction[idx_to_merge]

    # Target index is current index + attraction (e.g., 5 + (-1) = 4)
    target_idx <- idx_to_merge + att

    # 4. Merge
    # We remove the temp_id before merging to keep structure clean
    df$temp_id <- NULL
    df <- .merge_phases(df, idx_to_merge, target_idx)
  }

  # Final Cleanup
  # Sort by time
  df <- df[order(df$t_left), ]

  return(df)

}

##------------------------------------------------------------------------------------------------------------------
#' Merge phases
#' This function is heavily based on the implementation available on AMiGA.
#'
#' @param df A data.frame with the parameters of each potential growth phase.
#' @param idx1 row number one
#' @param idx2 row number two
#' @return A data frame with the merged growth phases.
#' @keywords internal
.merge_phases <- function(df, idx1, idx2){

  # 1. Identify rows to merge and sort them by time (t_left)
  # idx1 and idx2 are the row numbers in the current dataframe
  sub_df <- df[c(idx1, idx2), ]
  sub_df <- sub_df[order(sub_df$t_left), ]

  # 2. Compute aggregated values
  new_vals <- list(
    t_left     = min(sub_df$t_left, na.rm = TRUE),
    t_right    = max(sub_df$t_right, na.rm = TRUE),
    K          = sum(sub_df$K, na.rm = TRUE),
    r          = max(sub_df$r, na.rm = TRUE),
    r_left     = sub_df$r_left[1],         # The earliest left rate
    r_right    = tail(sub_df$r_right, 1),  # The latest right rate
    attraction = df$attraction[idx2]       # Keep attraction of the "target"
  )

  # 3. Update the DataFrame
  # We overwrite the second row (target) with new values
  # We assume columns match the order of new_vals list
  df[idx2, names(new_vals)] <- new_vals

  # 4. Remove the first row (source)
  df <- df[-idx1, ]

  # 5. Sort by time and reset row indices
  df <- df[order(df$t_left), ]
  rownames(df) <- NULL

  return(df)
}






##------------------------------------------------------------------------------------------------------------------
#' Pad inflection points and directions
#' This function is heavily based on the implementation available on AMiGA.
#'
#' @param ips A numeric vector of inflection points (indices).
#' @param its A numeric vector of inflection directions (1 for positive, -1 for negative).
#' @param edge A numeric value indicating which edge to pad (1 for left, -1 for right).
#' @param len A numeric value indicating the total length of the data vector (used for right edge padding).
#'
#' @return A list containing the padded inflection points and directions.
#'
#' @keywords internal
.pad_inflections <- function(ips, its, edge = 1, len = NULL) {

  # 1. Orient the data (Reverse if handling the Right Edge)
  # In R, we check if edge is -1
  if (edge == -1) {
    ips <- rev(ips)
    its <- rev(its)
    idx <- len  # In R, the last index is length (not length-1)
  } else {
    idx <- 1    # In R, the first index is 1 (not 0)
  }

  # 2. Logic Checks
  # Condition: If starts at index 1 AND is positive inflection
  if ((ips[1] == 1) && (its[1] == 1)) {
    # Do nothing (pass)

  } else if (its[1] == -1) {
    # If starts negative: Prepend [idx] and [1]
    ips <- c(idx, ips)
    its <- c(1, its)

  } else if (its[1] == 1) {
    # If starts positive (but not at index 1): Prepend [idx, idx] and [1, -1]
    ips <- c(idx, idx, ips)
    its <- c(1, -1, its)
  }

  # 3. Restore orientation (Reverse back if needed)
  if (edge == -1) {
    ips <- rev(ips)
    its <- rev(its)
  }

  # Return a named list
  return(list(ips = ips, its = its))
}


##------------------------------------------------------------------------------------------------------------------
#' Determine max growth
#'
#' @param df A data.frame with timepoint and mean fields.
#'
#' @return A data.frame with the maximum mean value and the timepoint where this occurs
#'
#'
#' @keywords internal
.determine_max_growth <- function(df){

  # Filter to the maximum mean
  max_growth <- df %>%

    # Determine max
    filter(mean == max(mean, na.rm = TRUE)) %>%

    # Format
    select(timepoint, mean) %>%

    rename("max_growth" = "mean",
           "max_growth_time" = "timepoint")


  return(max_growth)
}

##------------------------------------------------------------------------------------------------------------------
#' Determine auc
#'
#' @param df A data.frame with timepoint and mean fields.
#' @param tp_limit A numeric value indicating a timepoint t until which the AUC is calculated
#'
#' @return A data.frame with the empirical AUC from the mean GP fit
#'
#'
#' @importFrom pracma trapz
#' @keywords internal
.determine_auc <- function(df, tp_limit=NULL){

  if(is.null(tp_limit)){
    auc_value <- data.frame("AUC" = trapz(df$timepoint, df$mean))
  }else{

    # Rename column
    auc_colname = paste0("AUC_", tp_limit)

    # Filter until tp
    df_f <- df %>%
      filter(timepoint <= tp_limit)

    # AUC until tp
    auc_value <- data.frame("AUC" = trapz(df_f$timepoint, df_f$mean))
    colnames(auc_value) <- auc_colname

  }



  return(auc_value)
}


##------------------------------------------------------------------------------------------------------------------
#' Determine growth loss
#'
#' @param df A data.frame with timepoint and mean fields.
#' @param max_od_data A data.frame with the information output from .determine_max_growth()
#'
#' @return A data.frame with the total loss of growth from the moment when maximum growth is reached
#'
#' @keywords internal
.determine_growth_loss <- function(df, max_od_data){


  # Gather the data from the maximum onwards
  growth_loss_data <- df %>%

    # GAther data from the moment of max growth
    filter(timepoint >= max_od_data$max_growth_time) %>%
    mutate(growth_loss = mean - max_od_data$max_growth) %>%

    # Gather the most negative change
    filter(growth_loss == min(growth_loss, na.rm = TRUE)) %>%
    select(growth_loss)


  return(growth_loss_data)
}


#------------------------------------------------------------------------------------------------------------------
#' Estimate maximum growth rate
#'
#' @param df A data.frame with timepoint and first_derivative fields.
#' @return A data.frame with the maximum specific growth rate and doubling time
#'
#'
#' @keywords internal
.determine_max_growth_rate <- function(df, max_od_data){

 # Determine the maxium growth rate as the maximum of the first derivative
  max_gr <- df %>%
    #filter(timepoint <= max_od_data$max_growth_time) %>%
    filter(first_derivative == max(first_derivative, na.rm = TRUE)) %>%

    select(timepoint, first_derivative) %>%
    rename("max_growth_rate_timepoint" = "timepoint",
           "max_growth_rate" = "first_derivative") %>%

    mutate(doubling_time = log(2)/max_growth_rate)

  return(max_gr)
}


#------------------------------------------------------------------------------------------------------------------
#' Estimate maximum death rate
#'
#' @param df A data.frame with timepoint and first_derivative fields.
#' @param max_od_data A data.frame with the information output from .determine_max_growth()
#'
#' @return A data.frame with the maximum death rate
#'
#'
#' @keywords internal
.determine_max_death_rate <- function(df, max_od_data){



  # Determine the maximum death rate as the maximum of the first derivative
  min_gr <- df %>%
    # GAther data from the moment of max growth
    #filter(timepoint >= max_od_data$max_growth_time,
    #       !is.na(first_derivative))

    filter(!is.na(first_derivative))


  if(nrow(min_gr) == 0){

    min_gr <- data.frame("max_death_rate_timepoint" = NA,
                               "max_death_rate" = NA)


  }else{

    min_gr <- min_gr %>%
      filter(first_derivative == min(first_derivative, na.rm = TRUE)) %>%

      select(timepoint, first_derivative) %>%
      rename("max_death_rate_timepoint" = "timepoint",
             "max_death_rate" = "first_derivative")


  }

  return(min_gr)
}

#------------------------------------------------------------------------------------------------------------------
#' Estimate lag phase
#'
#' @param df A data.frame with timepoint and second_derivative fields.
#' @param max_gr_data A data.frame with the output of .determine_max_growth_rate().
#'
#' @return A data.frame with the timepoint at the end of lag-phase
#'
#' @keywords internal
.determine_lag_phase <- function(df, max_gr_data){

  # First only look at changes in rate before the timepoint where the maximum
  # growth rate occurs.
  lagphase_tp <- df %>%
    filter(timepoint <= max_gr_data$max_growth_rate_timepoint) %>%

    # Now find the point with the fastest acceleration. Likely the end of lag phase.
    filter(second_derivative == max(second_derivative, na.rm = TRUE)) %>%

    select(timepoint, mean) %>%
    rename("lag_phase_timepoint" = "timepoint",
           "growth_at_lag_phase" = "mean")

  return(lagphase_tp)
}



#------------------------------------------------------------------------------------------------------------------
#' Estimate stationary phase
#'
#' @param df A data.frame with timepoint and second_derivative fields.
#' @param max_gr_data A data.frame with the output of .determine_max_growth_rate().
#'
#' @return A data.frame with the timepoint at the end of lag-phase
#'
#' @keywords internal
.determine_stationary_phase <- function(df, max_gr_data, max_od_data){

  # First only look at changes in rate before the timepoint where the maximum
  # growth rate occurs.
  statphase_tp <- df %>%
    filter(timepoint <= max_od_data$max_growth_time) %>%
    filter(timepoint >= max_gr_data$max_growth_rate_timepoint,
           !is.na(second_derivative))

  # Check if there is even any data left for us to look at. Could be that the
  # curve never stops growing

  if(nrow(statphase_tp) == 0){

    statphase_tp <- data.frame("stationary_phase_timepoint" = NA,
                               "growth_at_stationary_phase" = NA)


  }else{

    # Now find the point with the fastest de-acceleration Likely the end start of staionary phase.
    statphase_tp <- statphase_tp %>%
      filter(second_derivative == min(second_derivative, na.rm = TRUE)) %>%

      select(timepoint, mean) %>%
      rename("stationary_phase_timepoint" = "timepoint",
             "growth_at_stationary_phase" = "mean")


  }


  return(statphase_tp)
}

#------------------------------------------------------------------------------------------------------------------
#' Calculate the inflection point and its tangent of a curve
#'
#' This function calculates the inflection point (the point of maximum slope)
#' and the tangent line at this point of a curve. The curve is defined by the
#' provided dataframe which contains timepoints and corresponding optical density
#' values.
#'
#' @param XX A dataframe containing the timepoints in a column named "timepoint"
#' and the corresponding optical density values in a column named "od".
#'
#' @return A list containing the inflection point and the tangent line data. The
#' inflection point is returned as a dataframe row from the input dataframe, with
#' an additional column "curve_id" labeled as "inflection_point". The tangent line
#' data is a dataframe with timepoints and corresponding optical density values on
#' the tangent line, with a column "curve_id" labeled as "inflection_tangent".
#'
#' @keywords internal
#' @export
.get_inflection_point <- function(XX) {

  stop("NOT IMPLEMENTED YET")
}





