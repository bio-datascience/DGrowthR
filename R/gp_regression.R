#------------------------------------------------------------------------------------------------------------------
#' Fit a Gaussian process model to data with laGP
#'
#' This function fits a Gaussian process model to the provided data. It is an internal function that
#' really only does the fit. Nothing else is done or returned by this function.
#'
#' @param X_data A data.frame containing the covariates used for model fitting.
#' @param y A vector containing the values to be predicted.
#' @param delete A logical value. If TRUE, the Gaussian process model will be deleted after fitting.
#'
#' @return The input object, with the Gaussian process information as a list.
#'
#' @keywords internal
.fit_lagp <- function(X, y, delete=FALSE) {
  
  # init lengthscale d and nugget g
  d <- laGP::darg(NULL, X)
  
  # init of nugget can be tricky for some curves
  g <- tryCatch({
    g <- laGP::garg(list(mle = TRUE), y)
    g
  },
  error=function(cond){
    message("DGrowthR is setting new nugget priors...")
    
    # centered and increased y
    r2mean <- mean(y-mean(y))^2
    
    # Check if r2mean is 0
    if(r2mean == 0){
      r2mean <- median(y-mean(y))^2
    }
    
    # If still 0 check max
    if(r2mean == 0){
      r2mean <- max(y-mean(y))^2
    }
    
    g <- list(mle=TRUE,
              
              # Start at 0.1 of variance of y or at machine eps
              start = max(0.1*var(y), sqrt(.Machine$double.eps)),
              
              # Start at 0.1 of variance of y or at machine eps
              min = min(0.1*var(y), sqrt(.Machine$double.eps)),
              max = var(y),
              
              # This is tricky as we don't have the params for this
              # laGP::garg sets a=3/2 and b as ~1276.323/mean((y-mean(y))^2)
              ab=c(3/2, 1276.323/r2mean))
    
    g
  })
  
  #print(g)
  # Start an anisotropic GP model
  gpsepi <- laGP::newGPsep(X, y, d = d$start, g = g$start, dK = TRUE)
  # Fit GP using MLE to estimate the best lengthscale and nugget
  mle <- laGP::mleGPsep(gpsepi, param = "both", tmin = c(d$min, g$min), tmax = c(d$max, g$max))
  # Compute the log-likelihood
  llik <- laGP::llikGPsep(gpsepi)
  # Compute posterior
  posterior <- laGP::llikGPsep(gpsepi, dab = d$ab, gab = g$ab)
  
  # Do not include gpsep obj by default
  fit_info <- list(X = X, y = y, mle = mle, llik = llik, posterior = posterior)
  
  if (delete == TRUE) {
    laGP::deleteGPsep(gpsepi)
    fit_info$gpsepi <- NULL
  } else if (delete == FALSE) {
    fit_info$gpsepi <- gpsepi
  }
  
  return(fit_info)
}


#------------------------------------------------------------------------------------------------------------------
#' Fit a Gaussian process model to data 
#'
#' This function fits a Gaussian process model to the provided data.
#'
#' @param X_data A data.frame containing the covariates used for model fitting. Should have at least timepoint, timepoint_n, and od columns.
#' @param delete A logical value. If TRUE, the Gaussian process model will be deleted after fitting.
#'
#' @return The input object, with the Gaussian process model information as a list
#'
#' @keywords internal
.fit_gp <- function(X_data, delete=FALSE) {
  
  
  # Prepare X matrix
  X_matrix <- X_data %>% 
    
    # In principle everything else should be timepoint and covariates
    select(-c(timepoint_n, od)) %>% 
    as.matrix()
  
  # Gather od data for modelling
  y_vector <- X_data$od
  
  # Fit GP model
  fitted_lagp <- .fit_lagp(X_matrix, y_vector, delete=FALSE)
  
  # Return fitted model
  return(fitted_lagp)
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
   predicted_lagp <- .predict_lagp(XX, lagp_fit, complete_sigma=complete_sigma, delete=delete)
   
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
       predicted_lagp$prediction$prediction_dataframe <- .estimate_derivatives(predicted_lagp$prediction$prediction_dataframe)
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
#' 
#' @return An object with the fitted model and prediction results. The results include the predictions (p) and the prediction inputs (XX).
#'
#' @rdname fit_predict_gp
#' @export
fit_predict_gp <- function(df, t_steps=1000, complete_sigma=FALSE, prepare_dataframe=FALSE, estimate_derivatives=FALSE, delete=TRUE, pred_for_t=NULL){
  
  # Check if curve_id field is present
  if("curve_id" %in% colnames(df)){
    df <- df %>% 
      select(-curve_id)
  }
  
  if("contrast_val" %in% colnames(df)){
    contrast_value <- unique(df$contrast_val)
  }
  
  # Fit the object
  gpobj <- .fit_gp(df, delete = FALSE)
  
  
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
      
      new_timepoints <- sort(c(pred_for_t, XX$timepoint))
      
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
