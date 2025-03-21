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
#'
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
                                                  n_cores=1) {
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
                        sample_n_curves=100, save_gp_data = FALSE, n_cores=1){
    
    
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
      object <- .fit_and_estimate(object,  model_covariate, predict_n_steps, downsample_every_n_timepoints, sample_posterior_gpfit, od_auc_at_t, sample_n_curves, save_gp_data, n_cores)
      
    # Check if paramters for GP fit have changed  
      }else if(object@gpfit_info$gpfit_parameters$model_covariate != model_covariate || object@gpfit_info$gpfit_parameters$predict_n_steps != predict_n_steps || object@gpfit_info$gpfit_parameters$downsample_every_n_timepoints != downsample_every_n_timepoints){
  
      message("[DGrowthR::estimate_growth_parameters] >> Fitting new GP models due to change in hyperparameters...")
      object <- .fit_and_estimate(object,  model_covariate, predict_n_steps, downsample_every_n_timepoints, sample_posterior_gpfit, od_auc_at_t, sample_n_curves, save_gp_data, n_cores)
  
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
                              sample_posterior_gpfit, od_auc_at_t, sample_n_curves, save_gp_data, n_cores){
  
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
  
  # Gather downsampled data
  #object.ds <- downsample_timepoints(object, downsample_every_n_timepoints)
  
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
  estimate_list <- foreach::foreach(i = model_groups, .combine = c, .multicombine = TRUE, .init = c(),
                                    .options.snow = list(progress = progress)
  ) %dopar% {
    
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
        
        
        #.fit_gp(od_for_gp, downsample_every_n_timepoints=downsample_every_n_timepoints, delete = FALSE)
        
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
      growth_params <- .characterize_growth(gpfit$prediction$prediction_dataframe, od_auc_at_t)
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
    
    slot(object, "growth_parameters") <- bind_rows(estimate_list[names(estimate_list) == "growth_parameters"]) 
    
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
#' 
#' @return A data.frame with the requested growth parameters for the curve given as input
#'
#'
#' @keywords internal
.characterize_growth <- function(df, od_auc_at_t=NULL){

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
  
  return(gparams)
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
    filter(timepoint <= max_od_data$max_growth_time) %>% 
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
    filter(timepoint >= max_od_data$max_growth_time,
           !is.na(first_derivative))
  
  
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





