
#' #------------------------------------------------------------------------------------------------------------------
#' @title Plot FPCAs
#'
#' @description This function plots the first two principal components calculated
#'
#' @param object  DGrowthR object.
#' @param color Covariate by which to color the points in the scatterplot. Should be one of the columns in `metadata`. Defaults to NULL
#'
#' @return A scatterplot of the first two FPCs
#'
#' @export
setGeneric("plot_fpca", function(object, color = NULL) {
  standardGeneric("plot_fpca")
})

#' @rdname plot_fpca
#' @export
setMethod(
  f = "plot_fpca",
  signature = "DGrowthR",
  definition = function(object, color = NULL) {
    # Make sure that fpca has been done
    if (length(object@fpca) == 0) stop(paste("First estimate functional PCA with estimate_fpca()"))

    # Gather the output for FPCA
    fpca_obj <- object@fpca$fdapace_obj
    metadata_df <- object@metadata

    # Prepare fpca dataframe
    fpcaX <- fpca_obj$xiEst
    colnames(fpcaX) <- paste0("fpc", 1:ncol(fpcaX))

    fpca_df <- data.frame(fpcaX)
    fpca_df$curve_id <- names(fpca_obj$inputData$Lt)

    # Prepare axis names
    fp1_name <- paste0("FPC-1 (", round(fpca_obj$cumFVE[1] * 100, 2), "% explained var.)")
    fp2_name <- paste0("FPC-2 (", round((fpca_obj$cumFVE[2] - fpca_obj$cumFVE[1]) * 100, 2), "% explained var.)")


    # Plot the first two components
    if (is.null(color)) {
      fp.g <- ggplot(fpca_df, aes(x = fpc1, y = fpc2)) +
        geom_point(color = "darkgrey", alpha = 0.5) +
        theme_light() +
        labs(
          x = fp1_name,
          y = fp2_name
        )
    } else if (!is.null(color)) {
      
      if(color == "cluster_assignment"){
        
        color_metadata <- object@cluster_assignment %>% 
          rename("color_covar" = "cluster")
        
      }else{
        
        # Get the relevant metadata
        color_metadata <- metadata_df %>%
          select(all_of(c("curve_id", color))) %>%
          rename("color_covar" = all_of(color))
        
      }

      # Add metadata
      fpca_df <- fpca_df %>%
        left_join(color_metadata, by = "curve_id")

      fp.g <- ggplot(fpca_df, aes(x = fpc1, y = fpc2, color = color_covar)) +
        geom_point(alpha = 0.5) +
        theme_light() +
        labs(
          x = fp1_name,
          y = fp2_name,
          color = color
        )
    }

    return(fp.g)
  }
)

#------------------------------------------------------------------------------------------------------------------
#' @title Plot UMAP components
#'
#' @description This function plots the first two UMAP components calculated
#'
#' @param object  Degrowth object.
#' @param color Covariate by which to color the points in the scatterplot. Should be one of the columns in `metadata`. Defaults to NULL
#'
#' @return A scatterplot of the first two UMAP components
#'
#' @export
setGeneric("plot_umap", function(object, color = NULL) {
  standardGeneric("plot_umap")
})

#' @rdname plot_umap
#' @export
setMethod(
  f = "plot_umap",
  signature = "DGrowthR",
  definition = function(object, color = NULL) {
    # Make sure that umap has been done
    if (is.null(object@umap_coord)) stop(paste("First estimate functional UMAP with estimate_umap()"))

    # Gather the output for FPCA
    umap_df <- object@umap_coord
    metadata_df <- object@metadata


    # Plot the first two components
    if (is.null(color)) {
      fp.g <- ggplot(umap_df, aes(x = V1, y = V2)) +
        geom_point(color = "darkgrey", alpha = 0.5) +
        theme_light() +
        labs(
          x = "UMAP 1",
          y = "UMAP 2"
        )
    } else if (!is.null(color)) {
      
      if(color == "cluster_assignment"){
        
        color_metadata <- object@cluster_assignment %>% 
          rename("color_covar" = "cluster")
      }else{
        
        # Get the relevant metadata
        color_metadata <- metadata_df %>%
          select(all_of(c("curve_id", color))) %>%
          rename("color_covar" = all_of(color))
        
      }

      # Add metadata
      umap_df <- umap_df %>%
        left_join(color_metadata, by = "curve_id")

      fp.g <- ggplot(umap_df, aes(x = V1, y = V2, color = color_covar)) +
        geom_point(alpha = 0.5) +
        theme_light() +
        labs(
          x = "UMAP 1",
          y = "UMAP 2",
          color = color
        )
    }

    return(fp.g)
  }
)


#------------------------------------------------------------------------------------------------------------------
#' @title Plot Optical Density
#'
#' @description This function plots the raw optical density curves, with the possibility of specifying color and facet variables
#'
#' @param object An object of class "DGrowthR"
#' @param color A column of metadata to color growth curves
#' @param facet A column of metadata to facet growth curves plots
#' @param data A character value indicating which data should be plotted. If "raw_od" then unprocessed data is plotted. If "od_data" (default) then the pre-processed data is plotted if present.
#'
#' @return A ggplot object showing lineplots of growth curves.
#'
#'
#' @export
setGeneric("plot_growth_curves", function(object, color=NULL, facet=NULL, data="od_data") {
  standardGeneric("plot_growth_curves")
})

#' @rdname plot_growth_curves
#' @export
setMethod(f = "plot_growth_curves",
          signature = "DGrowthR",
          definition = function(object, color=NULL, facet=NULL, data="od_data") {

            # Gather OD data
            if(data == "raw_od"){
              # GAther raw od
              od_data <- object@raw_od
            }else{
              od_data <- object@od_data 
            }
            
            # Determine if data has been log transformed
            if(object@log_od){
              y_axis_title <- "log(OD)"
            }else{
              y_axis_title <- "OD"
            }
            
            # Gather metadata
            metadata <- object@metadata

            # If color not null update OD data
            if(!is.null(color)){
            
              if(color=="cluster_assignment"){
                metadata_color <- object@cluster_assignment %>%
                  rename("color_var" = "cluster") %>%
                  select(curve_id, color_var)
              }else{
                
                # Change and selec color variable
                metadata_color <- metadata %>%
                  rename("color_var" = all_of(color)) %>%
                  select(curve_id, color_var)
                
              }

              # Join color variable to od_Data
              od_data <- od_data %>%
                left_join(metadata_color, by="curve_id")
            }

            # If facet is not null update OD data
            if(!is.null(facet)){
          
              if(facet == "cluster_assignment"){
                
                #Change and selec facet variable
                metadata_facet <- object@cluster_assignment %>%
                  rename("facet_var" = "cluster") %>%
                  select(curve_id, facet_var)
                
                
              }else{
                
                
                #Change and selec facet variable
                metadata_facet <- metadata %>%
                  rename("facet_var" = all_of(facet)) %>%
                  select(curve_id, facet_var)
                
                
              }

              # Join color variable to od_Data
              od_data <- od_data %>%
                left_join(metadata_facet, by="curve_id")

            }


            # At this point 'od_data' should have all that is necessary
            if(is.null(color)){

              p <- ggplot(od_data, aes(x=timepoint, y=od, group=curve_id)) +
                geom_line(color="#C5C5C5", alpha=0.25)

            }else{

              p <- ggplot(od_data, aes(x=timepoint, y=od, group=curve_id, color=color_var)) +
                geom_line(alpha=0.25) +
                labs(color=color)

            }


            # If there is a facet variable then facet
            if(!is.null(facet)){
              p <- p + facet_wrap(~facet_var)
            }


            # Update theme and return
            p <- p + labs(y=y_axis_title, x="Time") + theme_light()

            return(p)

          })


#------------------------------------------------------------------------------------------------------------------
#' @title Plot Gaussian Process
#'
#' @description This function plots a specific gaussian process with the 95% confidence interval.
#'
#' @param object An object of class "DGrowthR"
#' @param gpfit_id The name of the GP model that should be plotted. Should be in the gpfit_info$gpfit_data
#' @param show_input_od A logical value indicating if the od curves used to create the GP should be plotted
#'
#' @return A ggplot object showing lineplots of growth curves and GP mean.
#'
#'
#' @export
setGeneric("plot_gp", function(object, gpfit_id=NULL, show_input_od=FALSE) {
  standardGeneric("plot_gp")
})

#' @rdname plot_gp
#' @export
setMethod(f = "plot_gp",
          signature = "DGrowthR",
          definition = function(object, gpfit_id=NULL, show_input_od=FALSE) {
            
            gpfit_id_str <- gpfit_id
            
            
            
            # Check that there is an input
            if(is.null(gpfit_id_str)){
              stop("Please provide a valid gpfit_id")
            }
            
            # Gather the GPfit data
            gpfit_data <- object@gpfit_info$gpfit_data
            gpfit_metadata <- object@gpfit_info$gpfit_parameters$model_covariate
            
            # Check that the gpfit_id exists
            if(!gpfit_id_str %in% gpfit_data$gpfit_id){
              stop("The provided gpfit_id is not present in the current modelled data.")
            }
            
            # Gather the relevant gp data
            gp_data <- gpfit_data %>% 
              filter(gpfit_id_str==all_of(gpfit_id)) %>% 
              # We have information about the variance of each timepoint, use this to get the 95% interval
              mutate(q1 = qnorm(0.05, mean = mean, sd = sqrt(diag.sigma)),
                     q2 = qnorm(0.95, mean = mean, sd = sqrt(diag.sigma)))
            
            # Start a ggplot
            p <- ggplot()
            
            # Gather the od data if requested
            if(show_input_od){
              od_data <- gather_od_data(object, metadata_field = gpfit_metadata, field_value = gpfit_id_str)
              
              # Start by plotting the original data
              p <- p +
                geom_line(data=od_data, aes(x=timepoint, y=od, group=curve_id), color="grey", linewidth=0.5)
            }
            
            # Now plot the gpfit
            p <- p + 
              geom_line(data=gp_data, aes(x=timepoint, y=mean), linewidth=2, color="black") + 
              geom_ribbon(data=gp_data, aes(x=timepoint, y=mean, ymin = q1, ymax = q2), fill = "#C5C5C5", color="black", linetype = "longdash", alpha = 0.2) +
              
              theme_classic() +
              labs(x="Time", 
                   y="Growth",
                   title = paste("Gaussian Process fit for", gpfit_id_str),
                   subtitle = "Mean function and 95% confidence interval")
            
            
            return(p)
    })


#------------------------------------------------------------------------------------------------------------------
#' @title Plot Growth Comparison
#'
#' @description This function plots the null and the alternative fits to the model
#'
#' @param object An object of class "DGrowthR"
#' 
#' @return A list with two ggplot objects. One for the null model and another one for the alternative.
#'
#'
#' @export
setGeneric("plot_growth_comparison", function(object) {
  standardGeneric("plot_growth_comparison")
})
#' @rdname plot_growth_comparison
#' @export
setMethod(f = "plot_growth_comparison",
          signature = "DGrowthR",
          definition = function(object) {
            
            # Make sure a growth comparison has been performed
            if(length(object@growth_comparison) == 0){
              stop("[DGrowthR::plot_growth_comparison] >> First run a growth comparison with growth_comparison()")
            }
            
            # Gather the common growth data
            od_data <- object@growth_comparison$input_od
            
            # Start the baseline plot with the raw growth curves
            p <- ggplot() + 
              geom_line(data=od_data, aes(x=timepoint, y=od, color=covariate_value, group=curve_id), linewidth=0.5, alpha=0.3) +
              
              theme_light() +
              labs(x="Time", 
                   y="Growth",
                   title = paste("Gaussian Process fit"),
                   subtitle = "Mean function and 95% confidence interval")
            
            
            # Add information from the alternative model
            gpfit_data.alternative <- object@growth_comparison$gpfit_alternative %>% 
              mutate(q1 = qnorm(0.05, mean = mean, sd = sqrt(diag.sigma)),
                     q2 = qnorm(0.95, mean = mean, sd = sqrt(diag.sigma)))
            
            p.alt <- p +  
              geom_line(data=gpfit_data.alternative, aes(x=timepoint, y=mean, color=covariate_value, group=covariate_value), linewidth=2) +
              geom_ribbon(data=gpfit_data.alternative, aes(x=timepoint, y=mean, group=covariate_value, ymin = q1, ymax = q2, color=covariate_value), 
                          fill = "#C5C5C5", linetype = "longdash", alpha = 0.2)
            
            
            # Now for the null
            gpfit_data.null <- object@growth_comparison$gpfit_null %>% 
              mutate(q1 = qnorm(0.05, mean = mean, sd = sqrt(diag.sigma)),
                     q2 = qnorm(0.95, mean = mean, sd = sqrt(diag.sigma)))
            
            p.null <- p +
              geom_line(data=gpfit_data.null, aes(x=timepoint, y=mean), color="black", linewidth=2) +
              geom_ribbon(data=gpfit_data.null, aes(x=timepoint, y=mean, ymin = q1, ymax = q2),
                          color="black",
                          fill = "#C5C5C5", linetype = "longdash", alpha = 0.2)
            
            return(list("null" = p.null, "alternative" = p.alt))
})



#------------------------------------------------------------------------------------------------------------------
#' @title Plot Optical Density over Time
#' @description Calculates the optical density between the condition and control fitted curves for each time point, plots both models as well as the optical density difference between both fits.
#'
#' @param object An object of class "DGrowthR".
#' @return A ggplot object representing the line plot.
#'
#' @export
setGeneric(
  name = "plot_od_difference",
  def = function(object) {
    standardGeneric("plot_od_difference")
  }
)

#' @rdname plot_od_difference
#' @export
setMethod(
  f = "plot_od_difference",
  signature = "DGrowthR",
  definition = function(object) {
    # Generate predicted data
    data.predicted.alternative_1 <- alternative(object) %>% 
      filter(covariate_value == object@growth_comparison$comparison_info[2])
    data.predicted.alternative_0 <- alternative(object) %>% 
      filter(covariate_value == object@growth_comparison$comparison_info[3])
    
    # Select and rename columns
    data.predicted.alternative_1 <- data.predicted.alternative_1 %>%
      select(timepoint, mean) %>%
      rename(od_alt = mean)
    
    data.predicted.alternative_0 <- data.predicted.alternative_0 %>%
      select(timepoint, mean) %>%
      rename(od_null = mean)
    
    # Combine data and calculate od_diff
    data.combined <- full_join(data.predicted.alternative_1, data.predicted.alternative_0, by = "timepoint") %>%
      mutate(od_diff = abs(od_null - od_alt))
    
    # Reshape data to long format for plotting
    data.combined.long <- data.combined %>%
      pivot_longer(
        cols = c(od_alt, od_null, od_diff),
        names_to = "variable",
        values_to = "value"
      )
    
    plot <- ggplot(data.combined.long, aes(x = timepoint, y = value, color = variable)) +
      geom_line(size = 1) +
      labs(
        x = "Timepoint", y = "Growth",
        color = "Condition"
      ) +
      scale_color_manual(
        values = c(
          od_alt = "#F8766D",
          od_null = "#00BFC4",
          od_diff = "#00BA38"
        ),
        labels = c(
          od_alt = "Treatment Fit",
          od_null = "Control Fit",
          od_diff = "Absolute Difference"
        )
      ) +
      theme_light() +
      theme(legend.position = "bottom") +
      ggtitle("Difference of growth over Time")
    
    
    return(plot)
  }
)


#------------------------------------------------------------------------------------------------------------------
#' @title Plot Likelihood Ratios
#'
#' @description This function plots a histogram of permuted Likelihood Ratios with the actual Likelihood Ratio indicated by a vertical line.
#'
#' @param object An object of class "DGrowthR".
#' @param num_bins Number of bins to be used in the histogram.
#'
#' @return A ggplot object representing the histogram plot.
#' 
#' @export
setGeneric(name = "plot_likelihood_ratio", function(object, num_bins=50) standardGeneric("plot_likelihood_ratio"))

#' @rdname plot_likelihood_ratio
#' @export
setMethod(
  f = "plot_likelihood_ratio",
  signature = "DGrowthR",
  definition = function(object, num_bins) {
    lr_observed <- lrt_statistic(object)
    lr_permuted <- sort(perm_test_result(object), decreasing = TRUE)
  
    # df_observed_repeated <- data.frame(lr = lr_observed_repeated, type = "Actual Likelihood Ratio")
    df_lr_ob <- data.frame(lr = lr_observed, type = "Actual Likelihood Ratio")
    df_lr_pe <- data.frame(lr = lr_permuted, type = "Permuted Likelihood Ratio")
    
    df <- rbind(df_lr_ob, df_lr_pe)
    
    # Plot the histograms
    plot <- ggplot(df, aes(x = lr, fill = type)) +
      geom_histogram(bins = num_bins, alpha = 0.8, position = "identity") +
      geom_vline(aes(xintercept = lr_observed), color = "#F8766D", linetype = "solid", size = 1) +
      scale_fill_manual(values = c("Permuted Likelihood Ratio" = "#00BFC4", "Actual Likelihood Ratio" = "#F8766D")) +
      labs(
        x = "Likelihood Ratio", y = "Frequency",
        title = "Histogram of Permuted and Actual Likelihood Ratios",
        fill = "Type",
        color = "Line Type",
        linetype = "Line Type"
      ) +
      theme_light() +
      theme(legend.position = "bottom")
    
    return(plot)
  }
)