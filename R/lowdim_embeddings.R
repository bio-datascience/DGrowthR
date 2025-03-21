#' @title Functional UMAP
#'
#' @description This function estimates the UMAP coordinates using the uwot package
#'
#' @param object  DGrowthR object.
#' @param data A string indicating the source data to embed. At the moment only embeds "od_data"
#' @param ... Additional arguments to pass to the optns argument of uwot
#'
#' @return The result from umap function from uwot package
#' @seealso [uwot::umap()]
#'
#' @export
setGeneric("estimate_umap", function(object, data="od_data", ...) {
  standardGeneric("estimate_umap")
})
#' @rdname estimate_umap
#' @importFrom uwot umap
#' @export
setMethod(
  f = "estimate_umap",
  signature = "DGrowthR",
  definition = function(object, data="od_data", ...) {
    message(paste("[DGrowthR::estimate_umap]", ">>", "Estimating functional UMAP"))
    
    # Prepare data for FPCA
    if(data == "od_data"){
      message(paste("[DGrowthR::estimate_umap]", ">>", "Embedding optical density data"))
      
      od_data <- object@od_data
      od_data.wide <- od_data %>%
        pivot_wider(id_cols=curve_id, names_from=timepoint, values_from=od) %>%
        column_to_rownames("curve_id")
      
    }
    
    # UMAP
    od_umap <- umap(od_data.wide, init = "pca", verbose = F, ...)
    
    # Make dataframe
    umap_df <- as.data.frame(od_umap) %>%
      rownames_to_column("curve_id")
    
    
    # Return updated object
    slot(object, "umap_coord") <- umap_df
    
    return(object)
  }
)

#------------------------------------------------------------------------------------------------------------------
#' @title Functional Principal Component Analysis
#'
#' @description This function estimates the functional principal components of the complete dataset using the fdapace package
#'
#' @param object  DGrowthR object.
#' @param data A string indicating the source data to embed. At the moment only embeds "od_data"
#' @param ... Additional arguments to pass to the optns argument of FPCA
#'
#' @return The result from FPCA function from fdapace package
#'
#' @seealso [fdapace::FPCA()]
#'
#' @export
setGeneric("estimate_fpca", function(object, data="od_data", ...) {
  standardGeneric("estimate_fpca")
})

#' @rdname estimate_fpca
#' @export
setMethod(
  f = "estimate_fpca",
  signature = "DGrowthR",
  definition = function(object, data="od_data", ...) {
    message(paste("[DGrowthR::estimate_fpca]", ">>", "Estimating functional Principal Components"))
    
    
    if(data=="od_data"){
      message(paste("[DGrowthR::estimate_fpca]", ">>", "Embedding optical density data"))
      # Prepare data for FPCA
      od_data <- object@od_data
      fpca_input <- fdapace::MakeFPCAInputs(od_data$curve_id, od_data$timepoint, od_data$od)
      
    }
    
    # Estimate FPCA
    fpca_obj <- fdapace::FPCA(fpca_input$Ly, fpca_input$Lt, optns = list(...))
    
    # Return updated object
    slot(object, "fpca") <- list("fdapace_obj" = fpca_obj)
    
    return(object)
  }
)



#------------------------------------------------------------------------------------------------------------------

#' Determine optimal value of epsilon with elbow-detection
#'
#' This function determines the optimal value of epsilon for density-based clustering based on the sorted
#' K-nearest neighbor distance curve. The optimal epsilon is the 'elbow' of the curve and is found as the point
#' with the largest distance to the (imaginary) line defined by the first and last point.
#'
#' @param X numeric matrix of coordinates in UMAP space
#' @param k numeric value that indicates the K nearest neighoburs to be used
#'
#' @return numeric value corresponding to the optimal eps value
#'
#' @keywords internal
# internal method used for determining optimal eps

.determine_optimal_eps <- function(X, k) {
  
  # Estimate K-NN distance
  kn.df <- data.frame(kn_dist = sort(dbscan::kNNdist(x=X, k=k)), # Sorted K-NN distance
                      index = 1:nrow(X))
  
  # Determine reference points for imaginary line
  rp1 <- as.numeric(kn.df[1, ]) # Reference point 1 is the smallest K-NN distance
  rp2 <- as.numeric(kn.df[nrow(kn.df), ]) # Reference point 2 is the largest K-NN distance
  
  # Determine imaginary reference line difference
  rl <- rp2 - rp1
  
  # Determine distance of point to reference line
  kn.df$refline_dist <- apply(kn.df, 1, function(x){abs(det(cbind(rl, (x-rp1)))) /sqrt(sum(rl**2)) })
  
  
  # Gather point with largest distance to reference line
  optimal_eps <- kn.df %>%
    filter(refline_dist == max(refline_dist)) %>%
    
    # Gather optimal distance
    select(kn_dist) %>%
    unlist()
  
  
  # Return optimal eps
  return(optimal_eps)
}


#------------------------------------------------------------------------------------------------------------------
#' @title Clustering Function
#'
#' @description This function performs clustering using the specified algorithm (DBSCAN or GMM) on UMAP or fPCA coordinates.
#'
#' @param object A DGrowthR object containing preprocessed data including UMAP coordinates.
#' @param embedding A string indicating if the coordinates from "umap" or "fpca" should be used for clustering.
#' @param algorithm A string specifying the clustering algorithm to use ("dbscan" or "gmm").
#' @param k Numeric value indicating the number of nearest neighbors (for DBSCAN) or number of mixture components (for GMM).
#' If k < 1, it is interpreted as a proportion of the total number of curves.
#' @param eps Numeric value indicating the epsilon distance to use for DBSCAN clustering. If NULL, automatic selection based on k is performed.
#' @param ... Additional parameters to clustering algorithm
#'
#' @return Updated DGrowthR object with updated metadata dataframe indicating cluster membership.
#' @seealso [dbscan::dbscan()], [mclust::Mclust()]
#'
#' @export
setGeneric("clustering", function(object, ...) {
  standardGeneric("clustering")
})

#' @rdname clustering
#' @export
setMethod(
  f = "clustering",
  signature = "DGrowthR",
  definition = function(object, embedding="umap", algorithm = "dbscan", k = 0.01, eps = NULL) {
    
    # Determine embedding to be used
    if(embedding == "umap"){
      
      lowdim_coord <- object@umap_coord
      if(nrow(lowdim_coord) == 0) stop("First estimate functional UMAP with estimate_fpca()")
      message("[DGrowthR::clustering] >> Using UMAP coordinates")
    
    }else if(embedding == "fpca"){
      
      if(length(object@fpca) == 0) stop("First estimate functional UMAP with estimate_umap()")
      # Gather the output for FPCA
      fpca_obj <- object@fpca$fdapace_obj
      # Prepare fpca dataframe
      fpcaX <- fpca_obj$xiEst
      colnames(fpcaX) <- paste0("fpc", 1:ncol(fpcaX))
      
      lowdim_coord <- data.frame(fpcaX)
      lowdim_coord$curve_id <- names(fpca_obj$inputData$Lt)
      message("Clustering using fPCA coordinates...")
      
    }else{
      stop("Parameter 'embedding' should be one of 'fpca' or 'umap'.")
    }
    
   
    # Prepare lowdim data
    lowdim.data <- lowdim_coord %>%
      column_to_rownames("curve_id") %>%
      as.matrix()
    
    if (k < 1) {
      k <- round(k * nrow(lowdim.data))
    }
    
    if (algorithm == "dbscan") {
      message("[DGrowthR::clustering] >> Density based clustering with DBSCAN")
      message(paste("[DGrowthR::clustering] >> Using K-nearest neighbors:", k))
      if (is.null(eps)) {
        message("[DGrowthR::clustering] >> Automatic determination of optimal eps...")
        eps <- .determine_optimal_eps(lowdim.data, k)
      }
      message(paste("[DGrowthR::clustering] >> Using eps value:", eps))
      dbclusts <- dbscan::dbscan(lowdim.data, minPts = k, eps = eps)
      cluster.assignment <- data.frame(
        curve_id = rownames(lowdim.data),
        cluster = factor(dbclusts$cluster)
      )
      
    } else if (algorithm == "gmm") {
      message("[DGrowthR::clustering] >> Gaussian Mixture Model based clustering")
      message(paste("[DGrowthR::clustering] >> Using GMM with", k, "components"))
      gmm_result <- mclust::Mclust(lowdim.data, G = k)
      cluster.assignment <- data.frame(
        curve_id = rownames(lowdim.data),
        cluster = factor(gmm_result$classification)
      )
      
    } else {
      stop("Unsupported clustering algorithm. Please choose either 'dbscan' or 'gmm'.")
    }
    
    message(paste("[DGrowthR::clustering] >> Updating cluster_assignment slot.."))
    slot(object, "cluster_assignment") <- cluster.assignment
    
    message("[DGrowthR::clustering] >> Finished")
    return(object)
  }
)
