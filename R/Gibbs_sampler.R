## Gibbs sampling algorithm for 'SpaTopic'

theta_est<-function(Ndk, alpha){
  Nd<-apply(Ndk,1,sum)
  n.topics<-ncol(Ndk)
  apply(Ndk + alpha, 2, '/', Nd + n.topics * alpha)
}

beta_est<-function(Nwk, beta){
  Nk<-apply(Nwk,2,sum)
  n.words<-nrow(Nwk)
  apply(Nwk + beta, 1, '/', Nk + n.words * beta)
}

gibbs_spatial_LDA_multiple<-function(...){
  SpaTopic_inference(...)
}


#' 'SpaTopic': fast topic inference to identify tissue architecture in multiplexed images 
#' 
#' @description 
#' This is the main function of 'SpaTopic', implementing a Collapsed Gibbs
#' Sampling algorithm to learn topics, which referred to different tissue microenvironments, 
#' across multiple multiplexed tissue images. 
#' The function takes cell labels and coordinates on tissue images as input,
#' and returns the inferred topic labels for every cell, as well as topic contents, a distribution
#' over celltypes.
#' The function recovers spatial tissue architectures across images, 
#' as well as indicating cell-cell interactions in each domain.
#' 
#' @param tissue (Required). A data frame or a list of data frames. One for each image. 
#' Each row represent a cell with its image ID, X, Y coordinates on the image, celltype,
#' with column names (image, X, Y, type), respectively. For 3D tissue images, you may add 
#' either a 'Z' column (preferred) or 'Y2' column (legacy support) for the third dimension.
#' 
#' @param ntopics (Required). Number of topics. Topics will be obtained as distributions 
#' of cell types.
#' 
#' @param sigma Default is 50. The lengthscale of the Nearest-neighbor Exponential Kernel.
#' Sigma controls the strength of decay of correlation with distance in the kernel function.
#' Please check the paper for more information. 
#' Need to be adjusted based on the image resolution
#' 
#' @param region_radius Default is 400. The radius for each grid square when
#' sampling region centers for each image. 
#' Need to be adjusted based on the image resolution and pattern complexity.
#' 
#' @param npoints_selected Default is 1. Number of points sampled for each grid square 
#' when sampling region centers for each image. Used with \code{region_radius}.
#' 
#' @param kneigh Default is 5. Only consider the top 5 closest region centers for each cell.
#' 
#' @param ini_LDA Default is TRUE. Use warm start strategy for initialization and choose the best one
#' to continue. If 0, it simply uses the first initialization.
#' 
#' @param ninit Default is 10. Number of initialization. 
#' Only retain the initialization with the highest log likelihood (perplexity).
#' 
#' @param niter_init Default is 100. Warm start with 100 iterations in the Gibbs sampling 
#' during initialization.
#' 
#' @param beta Default is 0.05. A hyperparameter to control the sparsity of topic content
#' (topic-celltype) matrix \code{Beta}. A smaller value introduces more sparse in \code{Beta}.
#' 
#' @param alpha Default is 0.01. A hyperparameter to control the sparsity of document (region) content
#' (region-topic) matrix \code{Theta}. For our application, we keep it 
#' very small for the sparsity in \code{Theta}.
#' 
#' @param trace Default is FALSE. Compute and save log likelihood, \code{Ndk}, \code{Nwk} 
#' for every posterior samples. Useful when you want to use DIC to select number of 
#' topics, but it is time consuming to compute the likelihood for every posterior samples.
#' 
#' @param seed Default is 123. Random seed.
#' 
#' @param thin Default is 20. Key parameter in Gibbs sampling. 
#'  Collect a posterior sample for every thin=20 iterations.
#' 
#' @param burnin Default is 1000. Key parameter in Gibbs sampling.
#'  Start to collect posterior samples after 1000 iterations. You may increase
#'  the number of iterations for burn-in for highly complex tissue images.
#' 
#' @param niter Default is 200. Key parameter in Gibbs sampling. 
#' Number of posterior samples collected for model inference.
#' 
#' @param display_progress Default is TRUE. Display the progress bar. 
#' 
#' @param do.parallel Default is FALSE. Use parallel computing through R package \code{foreach}.
#' 
#' @param n.cores Default is 1. Number of cores used in parallel computing. 
#' 
#' @param axis Default is "2D". You may switch to "3D" for 3D tissue images. 
#' However, the model inference for 3D tissue is still under test. 
#' 
#' @param z_cellsize Default is region_radius*2. The thickness of each Z slice when
#' performing 3D stratified sampling. Only used when axis = "3D". Controls the 
#' Z-dimension binning resolution for region center selection in 3D tissue images.
#' Need to be adjusted based on the tissue thickness and Z-resolution.
#' 
#' @return Return a \code{\link{SpaTopic-class}} object. A list of outputs from Gibbs sampling. 
#' 
#' @seealso \code{\link{SpaTopic-class}}
#' 
#' @importFrom RANN nn2
#' 
#' @import foreach
#' 
#' @import methods
#' 
#' @import iterators
#' 
#' @examples 
#' 
#' ## tissue is a data frame containing cellular information from one image or
#' ## multiple data frames from multiple images.
#' 
#' data("lung5")
#' ## NOT RUN, it takes about 90s
#' library(sf)
#' #gibbs.res<-SpaTopic_inference(lung5, ntopics = 7,
#' #                               sigma = 50, region_radius = 400)
#'                              
#'                               
#' ## generate a fake image 2 and make an example for multiple images
#' ## NOT RUN
#' #lung6<-lung5
#' #lung6$image<-"image2"  ## The image ID of two images should be different
#' #gibbs.res<-SpaTopic_inference(list(A = lung5, B = lung6), 
#' #                 ntopics = 7, sigma = 50, region_radius = 400) 
#' 
#' @export
SpaTopic_inference<-function(tissue, ntopics, sigma = 50, region_radius = 400, kneigh = 5,
                                     npoints_selected = 1, ini_LDA = TRUE, ninit = 10, 
                                     niter_init = 100, beta = .05, alpha = .01,
                                     trace = FALSE, seed = 123, thin = 20, burnin = 1000,
                                     niter = 200, display_progress = TRUE, z_cellsize = region_radius*2,
                                     do.parallel = FALSE, n.cores = 1, axis = "2D"){
  
  set.seed(seed)
   
  if(is.data.frame(tissue)) tissue<-list(tissue)
  num_images<-length(tissue)
  
  if(num_images == 1) do.parallel<- 0 ## no need parallel computing if only one image
  
  ### combine multiple image into a single data frame
  itr_df<-do.call(rbind, tissue)
  
  ## check the colnames of the data frame
  if(!all(c("image","X","Y","type") %in% colnames(itr_df))){
    spatopic_message("ERROR", "Missing required columns. Please make sure you have image, X, Y, type in the colnames of tissue!")
    return(NULL)
  }
  
  if(any(is.na(itr_df$image)) | 
     any(is.na(itr_df$X)) |
     any(is.na(itr_df$Y)) |
     any(is.na(itr_df$type)) ){
    spatopic_message("ERROR", "Please make sure you have no NA in your input dataset!")
    return(NULL)
  }
  
  itr_df$X<-as.numeric(itr_df$X)
  itr_df$Y<-as.numeric(itr_df$Y)
  itr_df$type<-as.factor(itr_df$type) 
  itr_df$image<-as.factor(itr_df$image)
  
  if(length(levels(itr_df$image))<num_images){
    spatopic_message("ERROR", "Duplicate image ID! Please check images have distinct image ID!")
    return(NULL)
  }
  
  if(length(levels(itr_df$image))>num_images){
    spatopic_message("ERROR", "Multiple images should be input as a list of data frames 
                     or please check multiple image IDs in a single image")
    return(NULL)
  }
  
  
  if(axis == "3D" & ("Z" %in% colnames(itr_df) | "Y2" %in% colnames(itr_df))){
    if("Z" %in% colnames(itr_df)){
      itr_df$Z<-as.numeric(itr_df$Z)
    } else if("Y2" %in% colnames(itr_df)){
      itr_df$Y2<-as.numeric(itr_df$Y2)
    }
  }else{
    axis = "2D"
  }

  ## number of cells per image
  ncells<-table(itr_df$image)
  spatopic_message("INFO", paste("Number of cells per image:", paste(ncells, collapse = "\t")))
  
  
  ### coords for each sample
  coords<-lapply(tissue,GetCoords,axis = axis)
  rm(list= c("tissue"))

  spatopic_message("INFO", "Start initialization...")
  
  perplexity_min<-Inf
  Z_keep<-NULL
  D_keep<-NULL
  neigh_centers_keep<-NULL
  neigh_dists_keep<-NULL
  
  spatopic_message("INFO", paste("Number of Initializations:", ninit))
  
  ## if we could do parallel?
  is_doparallel_available <- requireNamespace("doParallel", quietly = TRUE)
  
  if(!is_doparallel_available & do.parallel){
    spatopic_message("WARNING", "Package 'doParallel' is not available. Running without parallel processing.")
    do.parallel <- FALSE
  }

  if(do.parallel){
    ### register cores for parallel computing 
    doParallel::registerDoParallel(n.cores)
  }
  
  for (ini in 1:ninit){ ### should chose the best one among random sample points
    
    #spatopic_message("INFO", paste("Initialization:", ini))

    i_id<-NULL  ## set global variable
    if(do.parallel){
      
      ### register cores for parallel computing 
      ## doParallel::registerDoParallel(n.cores)

      if(ini < 2){
        spatopic_message("INFO", paste("Parallel computing with number of cores:", n.cores))
      }
      
      results<-foreach::foreach(i_id = 1:num_images,.packages=c('sf','RANN'),
                       .export = c("stratified_sampling_sf", "stratified_sampling_3D_via_2D")) %dopar% {
        
        ## set seed for each core in parallel computing
        set.seed(seed+i_id)
        
        ## get coords
        coords_selected<-coords[[i_id]]
       
        if(axis == "3D" & "Z" %in% colnames(coords_selected)){
          ### 3D tissue image
          center_idx<-stratified_sampling_3D_via_2D(coords_selected, cellsize = c(region_radius*2,region_radius*2), z_cellsize = z_cellsize, num_samples_per_stratum = npoints_selected)
        }else{
          center_idx<-stratified_sampling_sf(coords_selected, cellsize = c(region_radius*2,region_radius*2), num_samples_per_stratum = npoints_selected)
        }
        
        ncenters<-length(center_idx)
        
        ### visualize selected centers
        coords_centers<-coords_selected[center_idx,]
        
        ### Find Knearest centers
        ##  Compute the kneigh nearest neighbors of each point in A with respect to B
        k_neighbor_centers <- nn2(coords_centers, coords_selected, k=kneigh, treetype = "kd",
                                  searchtype = "priority")
        neigh_centers<-k_neighbor_centers$nn.idx
        neigh_dists<-k_neighbor_centers$nn.dists
        
        return(list(neigh_centers = neigh_centers, neigh_dists = neigh_dists, 
                    coords_centers = coords_centers))
      }
    }else{
      
      results<-foreach(i_id = 1:num_images,.packages=c('RANN','sf')) %do% {
        
        ## get coords
        coords_selected<-coords[[i_id]]
    
        if(axis == "3D" & "Z" %in% colnames(coords_selected)){
          center_idx<-stratified_sampling_3D_via_2D(coords_selected, cellsize = c(region_radius*2,region_radius*2), z_cellsize = z_cellsize, num_samples_per_stratum = npoints_selected)
        }else{
          center_idx<-stratified_sampling_sf(coords_selected, cellsize = c(region_radius*2,region_radius*2), num_samples_per_stratum = npoints_selected)
        }
        
        ncenters<-length(center_idx)
        #print("number of centers selected:\n")
        #print(ncenters)
        
        ### visualize selected centers 
        ### [TODO] SAVE THE SELECTED CENTERS FOR EACH IMAGE
        ### [TODO] ADD OPTION TO MANUALLY SELECT CENTERS
        coords_centers<-coords_selected[center_idx,]
        
        ### Find Knearest centers
        ##  Compute the kneigh nearest neighbors of each point in A with respect to B
        k_neighbor_centers <- nn2(coords_centers, coords_selected, k=kneigh, treetype = "kd",
                                  searchtype = "priority")
        neigh_centers<-k_neighbor_centers$nn.idx
        neigh_dists<-k_neighbor_centers$nn.dists
        
        return(list(neigh_centers = neigh_centers, neigh_dists = neigh_dists, 
                    coords_centers = coords_centers))
      }
    }
  
  ## number of selected centers per image
  ncenters<-unlist(lapply(results, function(x) nrow(x$coords_centers)))
  #print("number of centers selected:")
  #print(ncenters)
  
  ## The nearest neighbor centers of each cells
  neigh_centers<-lapply(results, function(x) x$neigh_centers)
  neigh_centers<-do.call(rbind,neigh_centers)
  
  ### relabel centers if multiple images
  if(num_images > 1){
    add<-c(0,unlist(lapply(1:(length(ncenters)-1),function(x) sum(ncenters[1:x]))))
    neigh_centers<-neigh_centers+rep(add,times = ncells)
  }
    
  ## Distance of cells to their closest centers 
  neigh_dists<-lapply(results, function(x) x$neigh_dists)
  neigh_dists<-do.call(rbind,neigh_dists)
  
  ############ randomly sample Z and D ---------------------------------------------------
  
  ## Initialize Z, D for spatial Topic model
  neigh_centers <- neigh_centers - 1L
  D<-as.integer(neigh_centers[,1])  ### The closest region centers
  C<-as.integer(itr_df$type)-1L
  
  ### other parameters 
  ## number of topics
  K<-as.integer(ntopics)
  ## number of celltypes
  V<-length(unique((C)))
  ## number of regions
  M<-max(D)+1
  ## number of words
  N<-length(C)
  
  ######-----------------------------------------------------------------------------
  ## Both D and Z are changes over time
  
  ### additional input
  ## list of images and cells  
  docs<-as.matrix(cbind(as.integer(itr_df$image),C))
  ## list of different docs/regions
  doc_list<-as.integer(1:M)-1L
  ## list of different words
  word_list<-as.integer(1:V)-1L
  
    Z <- sample(1:ntopics,replace = TRUE, size = length(C))
    Z <- Z-1L
    
    Ndk <- table_2d_fast(D, Z, M, K) ## number of cells per topic per region (image specific)
    Nwk <- table_2d_fast(C, Z, V, K) ## number of cells per topic per celltype
    Nk <- table_1d_fast(Z, K)  ## number of cells in each topic 
    Nd <- table_1d_fast(D, M) ## number of cells in each region (image specific)
   
    
    ### initialization (warm start)
    if(ini_LDA){
    
      system.time(gibbs.res <-gibbs_sampler_c(docs, Ndk, Nwk, Nk, Nd, Z, D, neigh_dists, 
                                          neigh_centers, doc_list, word_list, 
                                          K, beta, alpha, sigma, thin, niter_init, 0,
                                          0, 0))
      perplexity<-gibbs.res$Perplexity
      
    }else{
      perplexity<-0
    }
    
    if(perplexity_min > perplexity){
      Z_keep<-Z
      perplexity_min <- perplexity
      D_keep<-D
      neigh_centers_keep<-neigh_centers
      neigh_dists_keep<-neigh_dists
    }
    ## within each iteration, release memory for the untouched objects
    ##gc()
  }
  
  ###  release memory for large items
  celltype<-levels(itr_df$type)
  cellname<-rownames(itr_df)
  rm(list= c("itr_df"))
  gc()
  
  ### Keep the best results during initialization
  if(ini_LDA){
    spatopic_message("RESULT", paste("Min perplexity during initialization:", round(perplexity_min, 4)))
  }
  Z<-Z_keep
  D<-D_keep
  neigh_centers<-neigh_centers_keep
  neigh_dists<-neigh_dists_keep

  reg_info <- paste("Number of region centers selected:", paste(ncenters, collapse = ", "))
  spatopic_message("INFO", reg_info)
  

  spatopic_message("INFO", paste("Average number of cells per region:", round(mean(table(D)), 2)))
  
  
  if(ninit > 1){ ## need update only when several initialization
     
      ## number of regions
    M<-max(D)+1  ## M might be changed since D changed 
    # A potential bug (what if the max(D) has been changed after gibbs sampling)
    
    Ndk <- table_2d_fast(D, Z, M, K) ## number of cells per topic per region (image specific)
    Nwk <- table_2d_fast(C, Z, V, K) ## number of cells per topic per celltype
    Nk <- table_1d_fast(Z, K)  ## number of cells in each topic 
    Nd <- table_1d_fast(D, M) ## number of cells in each region (image specific)
    
    ## list of different docs/regions
    doc_list<-as.integer(1:M)-1L
    ## list of different words
    word_list<-as.integer(1:V)-1L
  }
  
  spatopic_message("PROGRESS", "Initialization complete. Starting Gibbs sampling...")

  if(do.parallel){
    doParallel::stopImplicitCluster()
  }
  ## Gibbs_sampling 
  ## [TODO] should allow multiple chains to be run in parallel.
  ## Simply add results from multiple chains
  ## [TODO] may use Rcppparallel
  ## [TODO] compute log likelihood on a smaller subset of cells
  gibbs.res <-gibbs_sampler_c(docs, Ndk, Nwk, Nk, Nd, Z, D, neigh_dists, 
                                          neigh_centers, doc_list, word_list, 
                                          K, beta, alpha, sigma, thin, burnin, niter,
                                          trace, display_progress)
  
  gc()
  
  ## [TODO] make the output better and readable
  if(!trace){
    gibbs.res$doc.trace<-NULL
    gibbs.res$word.trace<-NULL
    gibbs.res$loglike.trace<-NULL
  }
  spatopic_message("COMPLETE", "Gibbs sampling completed successfully")
  spatopic_message("RESULT", paste("Final model perplexity:", round(gibbs.res$Perplexity, 4)))

   ## Model parameters used
  gibbs.res$parameters <- list(
    ## K
    beta = beta,
    alpha = alpha,
    sigma = sigma,
    region_radius = region_radius,
    kneigh = kneigh,
    npoints_selected = npoints_selected,
    thin = thin,
    burnin = burnin,
    niter = niter,
    ninit = ninit,
    ini_LDA = ini_LDA,
    niter_init = niter_init,
    trace = trace,            
    seed = seed,              
    display_progress = display_progress,  
    do.parallel = do.parallel,  
    n.cores = n.cores,          
    axis = axis,
    z_cellsize = z_cellsize
  )

  ## DIC  
  if(trace){
    Deviance<- -2*gibbs.res$loglike.trace
    gibbs.res$DIC<-0.5*var(Deviance)+mean(Deviance)
  }else{
    gibbs.res$DIC<-NULL
  }
  
  ## Beta: Topic Content
  gibbs.res$ntopics<-K
  gibbs.res$Beta<-as.data.frame(gibbs.res$Beta)
  colnames(gibbs.res$Beta)<-paste0("topic",1:K)
  rownames(gibbs.res$Beta)<-celltype
  
  ## Z.trace: posterior prob for every single cell
  gibbs.res$Z.trace<-as.data.frame(gibbs.res$Z.trace/niter)
  colnames(gibbs.res$Z.trace)<-paste0("topic",1:K)
  rownames(gibbs.res$Z.trace)<-cellname

  ## Cell topic assignments - most likely topic for each cell
  gibbs.res$cell_topics <- apply(as.matrix(gibbs.res$Z.trace), 1, which.max)

  ## Set class for the output
  class(gibbs.res) <- "SpaTopic"
  
  return(gibbs.res)
}

#' Print method for SpaTopic objects
#'
#' @description 
#' Provides a formatted summary of SpaTopic results when the object is printed.
#' This method displays key model metrics and explains how to access different components
#' of the model output.
#'
#' @param x An object of class "SpaTopic" returned by the SpaTopic_inference function
#' @param ... Additional arguments passed to print methods (not used)
#'
#' @details
#' The method displays:
#'   - Number of topics identified
#'   - Model perplexity (lower is better)
#'   - DIC (Deviance Information Criterion) for model comparison
#'   - A preview of the topic distributions across cell types
#'   - Instructions on how to access full results
#'
#' @return No return value, called for side effect of printing
#'
#' @examples
#' # If gibbs.res is a SpaTopic object:
#' # print(gibbs.res)
#'@export 
print.SpaTopic <- function(x, ...) {
  cat("SpaTopic Results\n")
  cat("----------------\n")
  # Display basic model information
  cat("Number of topics:", x$ntopics, "\n")
  cat("Perplexity:", x$Perplexity, "\n\n")  # Show model fit - lower is better
  # Only print DIC if it exists (which would be when trace=TRUE)
  if(!is.null(x$DIC)) {
    cat("DIC:", x$DIC, "\n\n")  # Deviance Information Criterion for model comparison
  }
  
  # Show preview of topic content
  cat("Topic Content(Topic distribution across cell types):\n")
  print(head(x$Beta, 5))
  cat("...\n\n")
  
  # Guide for accessing full results
  cat("Use $Z.trace for posterior probabilities of topic assignments for each cell\n")
  cat("Use $cell_topics for final topic assignments for each cell\n")
  cat("Use $parameters for accessing model parameters\n")
}

#' Format messages for SpaTopic package
#'
#' @description 
#' Creates consistently formatted messages for the SpaTopic package with 
#' timestamps and message type indicators. This function helps standardize
#' all output messages across the package. Error messages will stop execution.
#'
#' @param type Character string indicating message type (e.g., "INFO", "WARNING", "ERROR", "PROGRESS")
#' @param message The message content to display
#' @param timestamp Logical; whether to include a timestamp in the message (default: TRUE)
#'
#' @details
#' This function prefixes messages with a timestamp and the SpaTopic tag,
#' creating a consistent message format throughout the package.
#' When type="ERROR", this will stop execution with stop().
#' When type="WARNING", this will use warning() for non-fatal warnings.
#' All other message types will use message() for informational output.
#'
#' @return No return value, called for side effect of displaying a message
#'
#' @examples
#' \dontrun{
#' spatopic_message("INFO", "Starting analysis...")
#' spatopic_message("WARNING", "Parameter out of recommended range", timestamp = FALSE)
#' spatopic_message("ERROR", "Required input missing") # This will stop execution
#' spatopic_message("PROGRESS", "Processing complete")
#' }
#'
spatopic_message <- function(type = "INFO", message, timestamp = TRUE) {
  prefix <- paste0("[SpaTopic ", type, "]")
  if(timestamp) {
    time_str <- format(Sys.time(), "[%H:%M:%S]")
    prefix <- paste(time_str, prefix)
  }
  
  # Use appropriate function based on message type
  if(type == "ERROR") {
    stop(prefix, " ", message, call. = FALSE)
  } else if(type == "WARNING") {
    warning(prefix, " ", message, call. = FALSE)
  } else {
    message(prefix, " ", message)
  }
}

