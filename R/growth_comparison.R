#------------------------------------------------------------------------------------------------------------------
#' Predict the values of the fitted model to an object.
#'
#' This function fits a Gaussian Process model to a DGrowthR object and then generates predictions from the fitted model.
#' The function supports "alternative" and "null" models. The "alternative" model includes both time and contrast as predictors,
#' while the "null" model only includes time as a predictor.
#'
#' @param object A DGrowthR object to which the model should be fitted and then predictions should be generated.
#' @param comparison_info character vector describing the contrast to be made. The first value indicates the column in metadata where the contrast takes place. The third value is taken as the reference condition.
#' @param predict_n_steps A numeric value indicating the number of timepoints to make a prediction for with the fitted GP model. More points increases estimate accuracy, but adds compute time.
#' @param save_gp_data A logical indicating whether the fitted GP models should be stored for future analysis.
#' @param downsample_every_n_timepoints A numeric value indicating that the OD from every n timepoint should be used for GP fit. Might speed up fitting.
#' @param permutation_test A logical value indicating if a permutation test should be performed.
#' @param n_permutations A numerical values indicating the number of permutations to build in order to gather the null distribution of test statistics.
#' @param n_cores A numeric value indicating the number cores to be used. Useful when performing a permutation test.
#'
#' @return An object with the updated growth_comparison slot.
#'
#' @export
setGeneric(name = "growth_comparison", function(object,
                                                comparison_info,
                                                predict_n_steps=100,
                                                downsample_every_n_timepoints=1,
                                                save_gp_data=FALSE,
                                                permutation_test=FALSE,
                                                n_permutations=NULL,
                                                n_cores=1) standardGeneric("growth_comparison"))
#' @rdname growth_comparison
#' @export
setMethod(
  f = "growth_comparison",
  signature = "DGrowthR",
  definition = function(object,
                        comparison_info,
                        predict_n_steps=100,
                        downsample_every_n_timepoints=1,
                        save_gp_data=FALSE,
                        permutation_test=FALSE,
                        n_permutations=NULL,
                        n_cores=1) {

    message(paste("[DGrowthR::growth_comparison] >> Comparing", comparison_info[2], "to", comparison_info[3], "from the", comparison_info[1], "field."))

    # Gather downsampled data

    # Gather the necessary OD data
    od_data.comparison <- .gather_contrast(object, comparison_info, downsample_every_n_timepoints)


    # Fit the alternative GP
    gpfit.alternative <- fit_predict_gp(od_data.comparison,
                                        prepare_dataframe = TRUE,
                                        estimate_derivatives = FALSE,
                                        complete_sigma = FALSE,
                                        t_steps = predict_n_steps,
                                        delete=TRUE)

    # Prepare dataframe for null fit
    od_data.null <- od_data.comparison %>%
      select(-contrast_val)

    gpfit.null <- fit_predict_gp(od_data.null,
                                        prepare_dataframe = TRUE,
                                        estimate_derivatives = FALSE,
                                        complete_sigma = FALSE,
                                        t_steps = predict_n_steps,
                                        delete=TRUE)


    # Prepare output data.frame
    auc.treatment = .determine_auc(gpfit.alternative$prediction$prediction_dataframe %>% filter(contrast_val==1))$AUC
    auc.reference = .determine_auc(gpfit.alternative$prediction$prediction_dataframe %>% filter(contrast_val==0))$AUC

    max.treatment = .determine_max_growth(gpfit.alternative$prediction$prediction_dataframe %>% filter(contrast_val==1))$max_growth
    max.reference = .determine_max_growth(gpfit.alternative$prediction$prediction_dataframe %>% filter(contrast_val==0))$max_growth

    e.dist <- gpfit.alternative$prediction$prediction_dataframe %>%
      pivot_wider(id_cols = contrast_val, names_from = timepoint, values_from = mean) %>%
      column_to_rownames("contrast_val") %>%
      as.matrix() %>%
      dist() %>%
      as.vector()


    comparison_list <- list("result" = data.frame("comparison" = paste0(comparison_info[1], ": ", comparison_info[2], " v.s. ", comparison_info[3]),
                                    "likelihood_ratio" = gpfit.alternative$llik - gpfit.null$llik,
                                    "llik.alternative_model" = gpfit.alternative$llik,
                                    "llik.null_model" = gpfit.null$llik,

                                    "AUC.treatment" = auc.treatment,
                                    "AUC.reference" = auc.reference,
                                    "AUC.FoldChange" = auc.treatment/auc.reference,

                                    "max_growth.treatment" = max.treatment,
                                    "max_growth.reference" = max.reference,
                                    "max_growth.FoldChange" = max.treatment/max.reference,

                                    "euclidean.distance" = e.dist),

                            "input_od" = od_data.comparison %>%
                              mutate(covariate_value = if_else(contrast_val==1, comparison_info[2], comparison_info[3])) %>%
                              select(-contrast_val),

                            "comparison_info" = comparison_info

    )

    # Update the object slot
    slot(object, "growth_comparison") = comparison_list

    # Prepare predicted data
    if(save_gp_data){
      object@growth_comparison$gpfit_alternative <- gpfit.alternative$prediction$prediction_dataframe %>%
        mutate(covariate_value = if_else(contrast_val==1, comparison_info[2], comparison_info[3])) %>%
        select(-contrast_val)

      object@growth_comparison$gpfit_null <- gpfit.null$prediction$prediction_dataframe
    }

    if(permutation_test){
      object <- perm_test(object, num_perms = n_permutations,
                          n.cores = n_cores,
                          predict_n_steps=predict_n_steps,
                          gp.delete = TRUE)

      object@growth_comparison$result$empirical_p.value <- calculate_emp_p_value(object)
    }






     message("[DGrowthR::growth_comparison] >> Finished!")
     return(object)
  }
)

#------------------------------------------------------------------------------------------------------------------
#' @title Perform a permutation test for the DGrowthR object
#'
#' @description This function performs a permutation test for the DGrowthR object. It shuffles
#' the "contrast" column in the OD data and calculates the difference in log-likelihoods
#' between the alternative and null models for each permutation. The permutation test
#' helps assess the significance of the difference in models.
#'
#' @param object A DGrowthR object which contains the Gaussian Process regression results.
#' @param num_perms The number of permutations to perform. Defaults to 1000.
#' @param predict_n_steps A numeric value indicating the number of timepoints to make a prediction for with the fitted GP model. More points increases estimate accuracy, but adds compute time.
#' @param n.cores The number of CPU cores to use for parallel processing. Defaults to 1 (no parallel processing).
#' @param gp.delete A logical value indicating whether to delete the GPsep identifiers after the permutation test. Defaults to TRUE.
#' @param ... Additional arguments to be passed to the GP modelling function.
#'
#' @return A modified DGrowthR object with the calculated log-likelihood differences for each permutation.
#'
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @export
setGeneric(
  name = "perm_test",
  def = function(object, num_perms = 1000, n.cores = 1, predict_n_steps=100, gp.delete = TRUE, ...) {
    standardGeneric("perm_test")
  }
)
#' @rdname perm_test
#' @export
setMethod(
  f = "perm_test",
  signature = "DGrowthR",
  definition = function(object, num_perms = 1000, n.cores = 1, predict_n_steps=100, gp.delete = TRUE, ...) {

    message(paste("[DGrowthR::perm_test]", ">>", "Perform permutation test with", num_perms, "permutations"))


    # Parallel setup using forEach and dopar
    if (n.cores > 1) message(paste("[DGrowthR::perm_test]", ">>", "Register parallel backend with", n.cores, "cores"))
    # Source: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html
    cl <- snow::makeSOCKcluster(n.cores)
    doSNOW::registerDoSNOW(cl)

    od_perm <- object@growth_comparison$input_od
    pb <- txtProgressBar(max = num_perms, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)

    # To keep track of time
    time_started <- proc.time()[3]
    #doParallel::registerDoParallel(n.cores)

    # Perform permutations
    perms <- foreach::foreach(i = 1:num_perms, .combine = c, .multicombine = TRUE, .init = c(),
                              .options.snow = list(progress = progress)) %dopar% {
      # shuffle covariate label per curve_id
      perm_metadata <- od_perm %>%
        select(curve_id, covariate_value) %>%
        distinct() %>%
        mutate(covariate_value = sample(covariate_value, n(), replace=FALSE),
               contrast_val = if_else(covariate_value == object@growth_comparison$comparison_info[2], 1, 0)) %>%
        select(-covariate_value)


      # Gather the permuted growth data
      od_perm_shuffled <- od_perm %>%
        select(-covariate_value) %>%
        left_join(perm_metadata, by="curve_id")

      # Fit the alternative model
      gpfit.alternative.perm <- fit_predict_gp(od_perm_shuffled,
                                          prepare_dataframe = FALSE,
                                          estimate_derivatives = FALSE,
                                          complete_sigma = FALSE,
                                          delete=TRUE,
                                          t_steps = predict_n_steps)


      # Likelihood ratio from permutation test.
      return(gpfit.alternative.perm$llik)
    }

    close(pb)
    snow::stopCluster(cl)

    time_ended <- proc.time()[3]
    time_diff <- round(time_ended - time_started, 2)
    minutes <- floor(time_diff / 60)
    seconds <- round(time_diff %% 60)
    if (minutes > 0) {
      message(paste("[DGrowthR::perm_test] >> Time taken to perform permutation test:", minutes, "minute(s) and", seconds, "second(s)"))
    } else {
      message(paste("[DGrowthR::perm_test] >> Time taken to perform permutation test:", seconds, "second(s)"))
    }

    # Return the permuted ratio tests
    object@growth_comparison$permuted_llik.diff <- perms - object@growth_comparison$result$llik.null_model
    return(object)
  }
)

#------------------------------------------------------------------------------------------------------------------

#' Gather the data necesarry to perform a growth comparions.
#'
#' This function gathers the relevant data for the contrast indicated by the user. It is an internal function that
#' is used as part of the growth_comparison method for the DGrowthR class.
#'
#' @param object A DGrowthR object.
#' @param contrast_vector A 'vector' object containing three values: 1) The variable on which the contrast is made, 2) alternative condition, and 3) reference condition
#' @param downsample_every_n_timepoints A numeric value indicating that the OD from every n timepoint should be used for GP fit. Might speed up fitting.
#'
#' @return The input object, with the appropiate OD data.
#'
#' @keywords internal
# internal method used for gathering the relevant data
#.gather_contrast <- function(object, contrast_vector, downsample_every_n_timepoints) {

  # Gather the relevant data from contrast vector
#  metadata_field <- contrast_vector[1]
#  cov.alternative <- contrast_vector[2]
#  cov.reference <- contrast_vector[3]

  # Gather the relevant OD data
#  od_data.reference <- gather_od_data(object, metadata_field = metadata_field, field_value = cov.reference, downsample_every_n_timepoints=downsample_every_n_timepoints) %>%
#    mutate(contrast_val = 0)

#  od_data.alternative <- gather_od_data(object, metadata_field = metadata_field, field_value = cov.alternative, downsample_every_n_timepoints=downsample_every_n_timepoints) %>%
#    mutate(contrast_val = 1)




 # return(bind_rows(od_data.reference, od_data.alternative))
#}



#------------------------------------------------------------------------------------------------------------------

#' Compute the delta OD across all timepoints
#'
#' This function computes the delta || OD|| between two vectors
#'
#' @param vectorA A vector with one of the
#' @param contrast_vector A 'vector' object containing three values: 1) The variable on which the contrast is made, 2) alternative condition, and 3) reference condition
#' @param downsample_every_n_timepoints A numeric value indicating that the OD from every n timepoint should be used for GP fit. Might speed up fitting.
#'
#' @return The input object, with the appropiate OD data.
#'
#' @keywords internal
# internal method used for gathering the relevant data
.gather_contrast <- function(object, contrast_vector, downsample_every_n_timepoints=1) {

  # Gather the relevant data from contrast vector
  metadata_field <- contrast_vector[1]
  cov.alternative <- contrast_vector[2]
  cov.reference <- contrast_vector[3]

  # Gather the relevant OD data
  od_data.reference <- gather_od_data(object,
                                      metadata_field = metadata_field,
                                      field_value = cov.reference,
                                      downsample_every_n_timepoints = downsample_every_n_timepoints) %>%
    mutate(contrast_val = 0)

  od_data.alternative <- gather_od_data(object,
                                        metadata_field = metadata_field,
                                        field_value = cov.alternative,
                                        downsample_every_n_timepoints = downsample_every_n_timepoints) %>%
    mutate(contrast_val = 1)

  return(bind_rows(od_data.reference, od_data.alternative))
}
