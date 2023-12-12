## Gibbs sampling algorithm for SpaTopic

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


#' SpaTopic: fast topic inference to identify tissue architecture in multiplexed images 
#' 
#' @description 
#' This is the main function of SpaTopic, implementing a Collapsed Gibbs
#' Sampling algorithm to learn topics, which referred to different tissue microenvironments, 
#' across multiple multiplexed tissue images. 
#' The function takes as input cell labels and coordinates on tissue images 
#' and returns the inferred topic labels and topic contents, a distribution
#' over celltypes.
#' The function recovers spatial tissue architectures across images, 
#' as well as cell-cell interactions.
#' 
#' @param tissue (Required). A data frame or a list of data frames. One for each image. 
#' Each row represent a cell with its image ID, X, Y coordinates on the image, celltype,
#' with column names (image, X, Y, type), respectively. You may add another column 
#' Y2 for 3D tissue image.
#' 
#' @param ntopics (Required). Number of topics. Topics will be obtained as a distribution 
#' of 
#' 
#' @param sigma Default is 50. The lengthscale of the Nearest-neighbor Exponential Kernel.
#' Sigma controls the strength of decay of correlation with distance in the kernel function.
#' Please check the paper for more information. 
#' Need to be adjusted based on the image resolution
#' 
#' @param region_radius Default is 400. The radius for each grid square when
#' sampling region centers for each image. 
#' Need to be adjusted based on the resolution and the complexity in images.
#' 
#' @param npoints_selected Default is 1. Number of points sampled for each grid square 
#' when sampling region centers for each image. Used with \code{region_radius}.
#' 
#' @param kneigh Default is 5. Only consider the top 5 closest region centers for each cell.
#' 
#' @param ini_LDA Default is TRUE. Use warm start strategy for initialization and choose the best one
#' to continue. If 0, it simply just uses the first initialization.
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
#'  Collect a posterior sample for every 20 iterations.
#' 
#' @param burnin Default is 1000. Key parameter in Gibbs sampling.
#'  Start to collect posterior samples after 1000 iterations. You may increase
#'  the number of iterations for burn-in for highly complex tissue images.
#' 
#' @param niter Default is 200. Key parameter in Gibbs sampling. 
#' Number of posterior samples collected for the inference.
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
#' @return Return a \code{\link{gibbs.res-class}} object. A list of outputs from Gibbs sampling. 
#' 
#' @seealso \code{\link{gibbs.res-class}}
#' 
#' @importFrom RANN nn2
#' @import foreach
#' 
#' @examples 
#' 
#' ## tissue is a data frame containing cellular information from one image or
#' ## multiple data frames from multiple images.
#' 
#' data("lung5")
#' ## NOT RUN, it takes about 90s
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
                                     niter = 200, display_progress = TRUE,
                                     do.parallel = FALSE, n.cores = 1,axis = "2D"){
  
  set.seed(seed)
  
  if(is.data.frame(tissue)) tissue<-list(tissue)
  num_images<-length(tissue)
  
  if(num_images == 1) do.parallel<- 0 ## no need parallel computing if only one image
  
  ### combine multiple image into a single data frame
  itr_df<-do.call(rbind, tissue)
  
  ## check the colnames of the data frame
  if(!all(c("image","X","Y","type") %in% colnames(itr_df))){
    warning("Please make sure you have image, X, Y, type in the colnames of tissue!")
    return(NULL)
  }
  
  itr_df$X<-as.numeric(itr_df$X)
  itr_df$Y<-as.numeric(itr_df$Y)
  itr_df$type<-as.factor(itr_df$type) 
  itr_df$image<-as.factor(itr_df$image)
  
  if(length(levels(itr_df$image))<num_images){
    warning("Duplicate image ID! Please check images have distinct image ID!")
    return(NULL)
  }
  
  if(axis == "3D" & "Y2" %in% colnames(itr_df)){
    itr_df$Y2<-as.numeric(itr_df$Y2)
  }else{
    axis = "2D"
  }

  ## number of cells per image
  ncells<-table(itr_df$image)
  print("number of cells per image:")
  print(ncells)
  
  ### coords for each sample
  coords<-lapply(tissue,GetCoords,axis = axis)
  rm(list= c("tissue"))

  print("Initialization....")
  
  perplexity_min<-Inf
  Z_keep<-NULL
  D_keep<-NULL
  neigh_centers_keep<-NULL
  neigh_dists_keep<-NULL
  
  print("Numer of Initializations:")
  print(ninit)

  
  for (ini in 1:ninit){ ### should chose the best one among random sample points
    
    is_doparallel_available <- requireNamespace("doParallel", quietly = TRUE)
    
    if(!is_doparallel_available & do.parallel){
      warning("R package do.parallel is not available! The process is without Parallel.")
      do.parallel<-0
    }
    
    i_id<-NULL  ## set global variable
    if(do.parallel){
      
      ### register cores for parallel computing 
      doParallel::registerDoParallel(n.cores)
    
      results<-foreach::foreach(i_id = 1:num_images,.packages=c('sf','RANN'),
                       .export = c("stratified_sampling_idx_sf")) %dopar% {
        
        ## set seed for each core in parallel computing
        set.seed(seed+i_id)
        
        ## get coords
        coords_selected<-coords[[i_id]]
       
        ## Strategy 2 with library(sf)
        center_idx<-stratified_sampling_sf(coords_selected, 
                                               cellsize = c(region_radius*2,region_radius*2),
                                               npoints_selected)
        
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
    
        ## Strategy 2 with library(sf)
        center_idx<-stratified_sampling_sf(coords_selected, 
                                               cellsize = c(region_radius*2,region_radius*2),
                                               npoints_selected)
        
        ncenters<-length(center_idx)
        #print("number of centers selected:\n")
        #print(ncenters)
        
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
  M<-length(table(D))
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
  rm(list= c("itr_df"))
  gc()
  
  ### Keep the best results during initialization
  if(ini_LDA){
    print("Min perplexity during initialization:")
    print(perplexity_min)
  }
  Z<-Z_keep
  D<-D_keep
  neigh_centers<-neigh_centers_keep
  neigh_dists<-neigh_dists_keep

  print("number of region centers selected:")
  print(ncenters) 

  print("number of cells per region on average:")
  print(mean(table(D)))
  
  
  if(ninit > 1){ ## need update only when several initialization
     
      ## number of regions
    M<-length(table(D))  ## M might be changed since D changed
    
    Ndk <- table_2d_fast(D, Z, M, K) ## number of cells per topic per region (image specific)
    Nwk <- table_2d_fast(C, Z, V, K) ## number of cells per topic per celltype
    Nk <- table_1d_fast(Z, K)  ## number of cells in each topic 
    Nd <- table_1d_fast(D, M) ## number of cells in each region (image specific)
    
    ## list of different docs/regions
    doc_list<-as.integer(1:M)-1L
    ## list of different words
    word_list<-as.integer(1:V)-1L
  }
  
  ##print(proc.time())
  print("Finish initialization. Start Gibbs sampling....")
  
  ## Gibbs_sampling 
  ## [TODO] should allow multiple chains to be run in parallel.
  ## Simply add results from multiple chains
  ## [TODO] may use Rcppparallel
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
 
  print("Output model perplexity..")
  print(gibbs.res$Perplexity)
  
  return(gibbs.res)
}

