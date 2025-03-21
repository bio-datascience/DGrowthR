#' @title Preprocess Growth Data
#'
#' @description This function reads growth data from a file, applies several optional preprocessing steps, and saves the processed data.
#'
#' @param object A DGrowthR object containing the growth data.
#' @param baseline A numeric value specifying timepoint that will be used for baseline correction. If a vector is provided, then the minimum OD of the provided timepoints is used as baseline.
#' @param log_transform A logical value indicating whether to perform log transformation on optical density measurements. This enables accurate estimation of growth parameters.
#' @param skip_first_n_timepoints A numeric value the number of timepoints to skip for all growth curves. The analysis would start at skip_first_n_timepoints + 1.
#'
#' @return The processed DGrowthR object with updated od_data.
#'
#' @name preprocess_data
#' @rdname preprocess_data
#'
#' @export
setGeneric(name = "preprocess_data", function(object, baseline = 1, log_transform = TRUE, skip_first_n_timepoints = 0)
  standardGeneric("preprocess_data"))

#' @rdname preprocess_data
#' @export
setMethod(f = "preprocess_data",
          signature = "DGrowthR",
          definition = function(object, baseline = 1, log_transform = TRUE, skip_first_n_timepoints = 0){

            # Gather the raw data
            od_df_wide <- object@raw_od %>% 
              # Widen data 
              pivot_wider(id_cols = c(timepoint, timepoint_n), names_from = curve_id, values_from = od)

            # Assign unprocessed data to raw_od slot
            # slot(object, "raw_od") <- od_df

            # Validate inputs
            if (!is.numeric(baseline) && !is.null(baseline)) {
              stop("The 'baseline' argument must be a numeric vector specifying timepoint(s).")
            }

            if (!is.logical(log_transform)) {
              stop("The 'log_transform' argument must be a logical value.")
            }

            if (!is.numeric(skip_first_n_timepoints) && !is.null(skip_first_n_timepoints)) {
              stop("The 'skip_first_n_timepoints' argument must be a numeric value.")
            }
            
            # Pre-processing steps.
            
            ## Remove skipped timepoints
            if(skip_first_n_timepoints > 0){
              message(paste("Ignoring first", skip_first_n_timepoints, "for analysis."))
              
              od_df_wide <- od_df_wide %>% 
                # Filter timepoints to skip
                filter(timepoint_n > skip_first_n_timepoints)
            
            }
            
            # Re-pivot to make longer            
            od_df_long <- od_df_wide %>% 
              pivot_longer(cols = -c(timepoint, timepoint_n), names_to="curve_id", values_to="od") %>% 
              filter(!is.na(od))
            
            ## Log-transform data
            if(log_transform){
              
              # Check if any measurements are less than or equal to 0.
              if(any(od_df_long$od <= 0)){
                # If so, then gather a pseudocount for each growth curve equal to smallest non-zero od measurement
                message("Detected measurements of od <= 0. Adding pseudo od...")
                od_pseudocounts <- od_df_long %>% 
                  group_by(curve_id) %>% 
                  
                  # Only non-zero values
                  filter(od > 0) %>% 
                  
                  # Obtain the smallest non-zero OD
                  summarise(pseudo_od = min(od, na.rm = TRUE)) %>% 
                  ungroup()
                
                
                # Add the pseudo_od to each growth curve
                od_df_long %>% 
                  
                  left_join(od_pseudocounts, by="curve_id") %>% 
                  mutate(od = od + pseudo_od) %>% 
                  select(-pseudo_od)
              }
              
              # Perform log-transformation
              message("Log-transforming the growth measurements...")
              od_df_long <- od_df_long %>%
                mutate(od = log(od))
              
              
            }else if(!log_transform){
              message("Using linear growth measurements...")
              message("WARNING: Not log-transforming growth mesaurements might result in inaccurate estimation of growth parameters...")
            }
            
            
            ## Subtract baseline
            if(length(baseline) == 1){
              message(paste("Using OD from timepoint", baseline, "as baseline."))
            }else if(length(baseline)>1){
              str_timepoints <- paste0(baseline, collapse = ", ")
              message(paste("Using minimum OD from timepoints", str_timepoints, "as baseline"))
            }
            
            od_df_long <- od_df_long %>% 
              group_by(curve_id) %>% 
              mutate(od = od - min(od[timepoint_n %in% baseline]),
                     od = if_else(od < 0, 0, od)) %>% 
              ungroup()

            
            # Update OD data
            slot(object, "od_data") <- od_df_long
            
            # Update DGrowthR logged status
            slot(object, "log_od") <- log_transform
            
            # Update preprocessed status
            slot(object, "preprocessed") <- TRUE

            return(object)
          }
)