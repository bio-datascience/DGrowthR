#' @title Multiple Comparisons
#'
#' @description  This function processes multiple contrasts, performs an analysis, and returns a data frame with the results.
#'
#' Key features include:
#'
#' - **Multiple contrast processing**: Ability to process and analyze multiple contrasts.
#' - **Comprehensive Analysis**: Performs a range of analyses based on user preferences, including Likelihood Ratio, empirical p-value, area under the curve (AUC) for both condition and control, maximum growth rate of the condition, maximum difference, and maximum optical density (OD) for both condition and control.
#' - **Output and Display Formats**: Provides an option to save results to a TSV file.
#'
#' @param object DGrowthR object that has been instantiated.
#' @param comparison_list A list of of vectors. Each vector contains the information for an individual contrast.
#' @param predict_n_steps A numeric value specifying the number of time points for prediction. Defaults to 50.
#' @param downsample_every_n_timepoints A numeric value that indicates how timepoints are sampled to increase fit speed
#' @param n_permutations Number of permutations for the permutation test (default: 100).
#' @param permutation_test A logical value indicating if a permutation test should be performed.
#' @param n_cores Number of cores to use for the permutation test (default: 1).
#' @param write_to_tsv A logical value inndicating if the results should be save to a tsv file.
#' @param save_perm_stats A logical value indicating whether the permuted test statistics should be stored in the DGrowthR object.
#' @param filename Name of the TSV file.
#' @return Provides a data.frame summarizing the results.
#'
#'
#' @export
multiple_comparisons <- function(object, 
                                 comparison_list,
                                 predict_n_steps=50, 
                                 downsample_every_n_timepoints=1, 
                                 permutation_test=FALSE,
                                 n_permutations=NULL,
                                 n_cores=1,
                                 write_to_tsv=FALSE,
                                 save_perm_stats=FALSE,
                                 filename=NULL) {
  
  
  # Check if results should be written in TSV
  if(write_to_tsv){
    if(is.null(filename)){
      
      filename <- paste0("multiple_comparisons_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv")
      
    }
    
    message(paste("[DGrowthR::multiple_comparisons] >> Results will be written to", filename))
  }
  
  # Check that the contrast list is actually a list
  if (is.vector(comparison_list) && !is.list(comparison_list)) {
    comparison_list <- list(comparison_list)
  }
  
  # Iterate over contrasts and gather the results
  results_lists <- sapply(1:length(comparison_list),
                                
                                .multiple_comparisons, 
                                
                                object = object,
                                comparison_list=comparison_list,
                                predict_n_steps=predict_n_steps, 
                                downsample_every_n_timepoints=downsample_every_n_timepoints, 
                                save_gp_data=FALSE, 
                                permutation_test=permutation_test,
                                n_permutations=n_permutations,
                                n_cores=n_cores,
                                save_perm_stats = save_perm_stats,
                                write_to_tsv=write_to_tsv,
                                filename=filename)
  
  # Bind all test results
  if(save_perm_stats){
    result_df <- bind_rows(results_lists["result_df", ]) 
  } else{
    result_df <- bind_rows(results_lists)
  }
  
  # Adjust multiple tests
  result_df <- result_df %>% 
    mutate(pvalue.adjust = p.adjust(empirical_p.value, method = "BH"))
  
  if(write_to_tsv){
    write_tsv(result_df, filename)
  }
  
  
  # If return permutation
  if(save_perm_stats){
    
    perm_df <- bind_rows(results_lists["perm_values",])
    
    if(write_to_tsv){
      basefilename <- basename(filename)
      basepath <- dirname(filename)
      filenamep <- paste0(basepath, "/permutations_", basefilename)
      
      write_tsv(perm_df, filenamep)
    }
    
    return(list("result_df"=result_df,
                "permutation_df" = perm_df))
  }else{
    
    return(result_df)
    
  }
}



#------------------------------------------------------------------------------------------------------------------
#' A wrapper for performing multiple comparisons that can catch errors.
#'
#' @param object A DGrowthR object to which the model should be fitted and then predictions should be generated.
#' @param comparison_list A list of contrasts.
#' @param predict_n_steps A logical
#' @param downsample_every_n_timepoints A numeric value indicating that the OD from every n timepoint should be used for GP fit. Might seep up fitting.
#' @param save_gp_data A logical value, indicating if the mean GP values and GP fit parameters should be saved to object.
#' @param permutation_test A logical value indicating if a permutation test should be performed.
#' @param n_permutations A numerical values indicating the number of permutations to build in order to gather the null distribution of test statistics.
#' 
#' @keywords internal
.multiple_comparisons <- function(index,
                                  comparison_list,
                                  object,
                                  predict_n_steps=100, 
                                  downsample_every_n_timepoints=1, 
                                  save_gp_data=FALSE, 
                                  permutation_test=FALSE,
                                  n_permutations=NULL,
                                  n_cores=1,
                                  save_perm_stats=FALSE,
                                  write_to_tsv=FALSE,
                                  filename=NULL){
  
  # Gather the information for the contrast
  comparison_info <- comparison_list[[index]]
  contrast_str <- paste0(comparison_info[1], ": ", comparison_info[2], " vs ", comparison_info[3])
  message(paste("[DGrowthR::multiple_comparisons] >> Computing", contrast_str, "at step", index, "of", length(comparison_list)))
  
  result_list <- tryCatch({
    
    gobj <- growth_comparison(object, 
                      comparison_info = comparison_info,
                      predict_n_steps=predict_n_steps, 
                      downsample_every_n_timepoints=downsample_every_n_timepoints, 
                      save_gp_data=FALSE, 
                      permutation_test=permutation_test,
                      n_permutations=n_permutations,
                      n_cores=n_cores)
    
    if(save_perm_stats){
      
      # Prepare a dataframe
      perm_df <- data.frame(comparison = gobj@growth_comparison$result$comparison,
                            perm_vals = perm_test_result(gobj),
                            nperm = 1:n_permutations) %>% 
        pivot_wider(id_cols = comparison, names_from=nperm, values_from = perm_vals)
      
      
      list(result_df = gobj@growth_comparison$result,
           perm_values = perm_df)
      
    }else{
      list(result_df = gobj@growth_comparison$result)
    }
    
  },
  error = function(cond){
    message(paste("[DGrowthR::multiple_comparisons] >> Unable to perform comparison for:", paste0(comparison_info, collapse = ", ")))
    data.frame("comparison" = paste0(comparison_info[1], ": ", comparison_info[2], " v.s. ", comparison_info[3]))
  })
  
  if(write_to_tsv){
    write_tsv(result_list$result_df, filename, append = (index > 1))
    
  }
  
  return(result_list)
  
}

