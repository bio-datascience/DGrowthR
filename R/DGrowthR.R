#   _____   _____                   _   _     _____
#  |  __ \ / ____|                 | | | |   |  __ \
#  | |  | | |  __ _ __ _____      _| |_| |__ | |__) |
#  | |  | | | |_ | '__/ _ \ \ /\ / / __| '_ \|  _  /
#  | |__| | |__| | | | (_) \ V  V /| |_| | | | | \ \
#  |_____/ \_____|_|  \___/ \_/\_/  \__|_| |_|_|  \_\ by Medina Feldl & Roberto Olayo-Alarcon
#
# -------------------------------------------------------------------------------------------------------------
#' @title This S4 class represents a DGrowthR model.
#' @description It is used to store and validate data and parameters necessary for DGrowthR analysis.
#'
#' @slot od_data A data.frame expected to contain the following columns: "timepoint", "od", and "curve_id". These represent the time points, and optical density (od), and an identifier for each curve respectively.
#' @slot metadata A data.frame that contains the "curve_id" column along with any metadata variables
#' @slot raw_od Optical density information before pre-processing or re-formatting.
#' @slot verbose A logical value indicating whether the output should be suppressed. Defaults to TRUE, i.e., output is not suppressed.
#' @slot fpca An fpca object as the output of FPCA from fdapace
#' @slot umap_coord Coordinates for UMAP embedding
#' @slot cluster_assignment A data.frame matching curve_ids to the cluster assigned by clustering() function.
#' @slot gpfit_info A list with the estimated mean gp process and the parameters used for fitting.
#' @slot growth_comparison The results of a growth comparison.
#' @slot log_od A logical indicating if the OD data has been logged.
#' @slot preprocessed A logical variable indicating if the od data has been preprocessed with the preprocess_data() function.
#' @slot growth_parameters A data.frame with the estimated growth parameters after gp modelling
#'
#' @section Validity:
#' For an object to be a valid DGrowthR, it must meet the following criteria:
#' \itemize{
#'  \item {The "od_data" slot must be a data frame containing the columns "timepoint", "timepoint_n", "curve_id", and "od".}
#'  \item {The "metadata" slot must be a data frame with the "curve_id" column.}
#'  \item {All of the "curve_id"s in the "od_data" data frame must be present in the "curve_id" column of the metadata.}
#' }
#'
#' @export
DGrowthR <- setClass(
  "DGrowthR",
  slots = list(
    od_data = "data.frame",
    metadata = "data.frame",
    raw_od = "data.frame",
    verbose = "logical",
    fpca = "list",
    umap_coord = "data.frame",
    cluster_assignment = "data.frame",
    gpfit_info = "list",
    growth_comparison = "list",
    log_od = "logical",
    preprocessed = "logical",
    growth_parameters = "data.frame"
  ),
  
  prototype = list(od_data = data.frame(),
                   metadata = data.frame(),
                   raw_od = data.frame(),
                   verbose = FALSE,
                   fpca = list(),
                   umap_coord = data.frame(),
                   cluster_assignment = data.frame(),
                   gpfit_info = list(),
                   growth_comparison = list(),
                   log_od = FALSE,
                   preprocessed = FALSE,
                   growth_parameters = data.frame()),
  
  validity = function(object) {
    if (!("timepoint" %in% colnames(object@od_data))) {
      return("timepoint column not found in input data")
    }
    if (!("curve_id" %in% colnames(object@od_data))) {
      return("curve_id column not found in input data")
    }
    if (!("od" %in% colnames(object@od_data))) {
      return("od column not found in input data")
    }
    if (!("curve_id" %in% colnames(object@metadata))) {
      return("curve_id column not present in metadata")
    }
    if (!all(unique(object@od_data$curve_id) %in% unique(object@metadata$curve_id))) {
      return("Not all curve_id values in od_data are present in the metadata")
    }
    
    return(TRUE)
  }
  
)


#' @title Instiantiate a DGrowthR object.
#'
#' @description This function receives a dataframe with the optical density measurements, and, optionally and metadata  data frame.
#'
#' @param od_data A data.frame containing all of the optical density data. It should have a "Time" column with each row being a timepoint and the rest of the columns being an individual growth curve. The names of the growth curves should match those in the "curve_id" field in metadata.
#' @param metadata Optional. A data.frame containing a "curve_id" field and all additional metadata. If none is provided, then one is created automatically containing only the wells from the curve_ids in od_data.
#' @param verbose A logical variable indicating if DGrowthR object should be verbose.
#' 
#' @return A DGrowthR object with filled raw_od, od_data, and metadata slots.
#'
#' @name DGrowthRFromData
#' @rdname DGrowthRFromData
#' @export
DGrowthRFromData <- function(od_data, metadata=NULL, verbose=TRUE){
  
  # Check that a "Time" column exists in od_data
  if(!"Time" %in% colnames(od_data)){
    stop('No column called "Time" in od_data')
  }
  
  # Check that there are no repeated curve_ids
  if(ncol(od_data)-1 != length(unique(colnames(od_data)))-1){
    stop("There are repeated curve identifiers. Make sure all colnames are unique.")
  }
  
  # Check if metadata is provided, if not, then create a dumb metadata file
  if(is.null(metadata)){
    message("Creating metadata from od_data")
    
    # Gather all of the curve_ids
    curve_ids <- colnames(od_data)[colnames(od_data) != "Time"]
    metadata <- data.frame(curve_id = curve_ids,
                           well = str_split_i(curve_ids, "_", -1))
    
  }
  
  # Pivot od_data
  od_data_long <- od_data %>% 
    mutate(timepoint_n = 1:n()) %>% 
    pivot_longer(cols = -c(Time, timepoint_n), names_to = "curve_id", values_to = "od") %>% 
    rename("timepoint" = "Time")
  
  
  # Give some summary statistics,
  message(paste("Creating a DGrowthR object..."))
  message(paste(ncol(od_data)-1, "growth curves."))
  message(paste(nrow(od_data), "timepoints."))
  message(paste(ncol(metadata)-1, "covariates (including well)."))
  
  
  # Create the DGrowthR object
  object <- new("DGrowthR", 
                od_data = od_data_long, 
                metadata = metadata,
                raw_od = od_data_long,
                verbose=verbose)
  
  return(object)
  
}

#------------------------------------------------------------------------------------------------------------------
#' @title Gather specific growth data from the DGrowthR object
#' @description Extracts the requested growth curves according to the requested metadata fields and variables.
#'
#' @param object An object of class "DGrowthR".
#' @param metadata_field A character string specifying one of the fields in the metadata from which the variables will be searched in.
#' @param field_value The value of the metadata_field that for which the relevant growth curves will be extracted. 
#' @param downsample_every_n_timepoints A numeric value indicating that the OD from every n timepoint should be used for GP fit. Might seep up fitting.
#'
#' @return A data.frame containing the optical density data requested.
#'
#' @export
setGeneric(name = "gather_od_data", function(object, metadata_field, field_value, downsample_every_n_timepoints=1) standardGeneric("gather_od_data"))
#' @rdname gather_od_data
#' @export
setMethod(f = "gather_od_data", signature = "DGrowthR", definition = function(object, metadata_field, field_value, downsample_every_n_timepoints=1) {
  
  # Gather the od_data and metadata
  od_data <- object@od_data
  metadata <- object@metadata
  
  
  # If the requested field is curve_id, then just return that
  if(metadata_field == "curve_id"){
    
    req.data <- od_data %>% filter(curve_id == field_value)
    
  }else{
    # Gather the relevant growth curve ids
    curve_ids_requested <- metadata %>% 
      rename("metadata_field" = all_of(metadata_field)) %>% 
      
      filter(metadata_field == field_value) %>% 
      select(curve_id) %>% 
      unlist()
    
    req.data <- od_data %>% filter(curve_id %in% curve_ids_requested)
    
  }
  
  if (downsample_every_n_timepoints > 1) {
    timepoints_to_keep <- seq(from = min(req.data$timepoint_n), to = max(req.data$timepoint_n), by = downsample_every_n_timepoints)
    req.data <- req.data %>% filter(timepoint_n %in% timepoints_to_keep)
    
  }
  
  
  
  # Return the requested growth curves
  return(req.data)
  
  
})




# fix tidyr and dplyr errors during build
utils::globalVariables(c("Time", "V1", "V2", "color_covar", "color_var", "contrast_list", "contrast_val",
                         "covariate_value", "curve_id", "diag.sigma", "empirical_p.value", "facet_var",
                         "first_derivative", "fpc1", "fpc2", "gpfit_id", "growth_loss", "i", "kn_dist", "lr",
                         "max_growth_rate", "object.comparison", "od", "od_alt", "od_diff", "od_null", "pseudo_od",
                         "q1", "q2", "refline_dist", "sampled_id", "second_derivative", "str_remove",
                         "str_split_i", "timepoint", "timepoint_n", "type", "value", "variable", "well",
                         "perm_vals", "comparison", "nperm"
))


#' Accessor for the alternative slot
#'
#' This function retrieves the alternative slot from the provided DGrowthR object.
#'
#' @param object An object from which to extract the alternative slot.
#' @return The alternative slot of the object.
#' @export
setGeneric(name = "alternative", function(object) standardGeneric("alternative"))

#' @rdname alternative
#' @export
setMethod(f = "alternative", signature = "DGrowthR", definition = function(object) object@growth_comparison$gpfit_alternative)

#' Accessor for the null slot
#'
#' This function retrieves the null slot from the provided DGrowthR object.
#'
#' @param object An object from which to extract the null slot.
#' @return The null slot of the object.
#' @export
setGeneric(name = "null", function(object) standardGeneric("null"))

#' @rdname null
#' @export
setMethod(f = "null", signature = "DGrowthR", definition = function(object) object@growth_comparison$gpfit_null)

#' Accessor for the perm_test_result slot
#'
#' This function retrieves the perm_test_result slot from the provided DGrowthR object.
#'
#' @param object An object from which to extract the perm_test_result slot.
#' @return The perm_test_result slot of the object.
#' @export
setGeneric(name = "perm_test_result", function(object) standardGeneric("perm_test_result"))
#' @rdname perm_test_result
#' @export
setMethod(f = "perm_test_result", signature = "DGrowthR", definition = function(object) object@growth_comparison$permuted_llik.diff)


#' Accessor for the likelihood ratio test result
#'
#' This function retrieves the likelihood ratio test result from the growth comparison slot
#'
#' @param object An object from which to extract the perm_test_result slot.
#' @return The LRT test statistic 
#' @export
setGeneric(name = "lrt_statistic", function(object) standardGeneric("lrt_statistic"))
#' @rdname lrt_statistic
#' @export
setMethod(f = "lrt_statistic", signature = "DGrowthR", definition = function(object) object@growth_comparison$result$likelihood_ratio)


#' Accessor for the result of growth comparison
#'
#' This function retrieves the result of a growth comparison
#'
#' @param object An object from which to extract the perm_test_result slot.
#' @return The growth comparison result data.frame 
#' @export
setGeneric(name = "growth_comparison_result", function(object) standardGeneric("growth_comparison_result"))
#' @rdname growth_comparison_result
#' @export
setMethod(f = "growth_comparison_result", signature = "DGrowthR", definition = function(object) object@growth_comparison$result)




#------------------------------------------------------------------------------------------------------------------
#' Display a summary of the DGrowthR object
#'
#' This method displays a summary of the DGrowthR object, including the first few rows
#' of the od_data slot and the log-likelihood and posterior probabilities of the alternative
#' and null models.
#'
#' @param object A DGrowthR object.
#'
#' @export
setMethod(
  f = "show",
  signature = "DGrowthR",
  definition = function(object) {
    cat(is(object)[[1]], "\n", sep = "")
    cat("Head(Data):\n")
    print(head(object@od_data))
  }
)

