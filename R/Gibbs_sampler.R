#library(ggplot2)
#library(RANN)
#library(Rcpp)
#library(RcppArmadillo)
#library(RcppProgress)
#library(sp)
#library(foreach)
#library(scales)

#sourceCpp('src/Gibbs.cpp')
#source("R/spatial_func.R")

## Rewrite this in Rcpp
perplexity_spaLDA<-function(theta_dk, beta_wk, C, neigh_dists, sigma,
                            neigh_centers, n.doc, alpha, beta){
  
  n.cells<-length(C)
  Gaus<-exp(-neigh_dists/sigma)
  Gaus_sum<-apply(Gaus,1,sum)
  eta<-sweep(Gaus,1,Gaus_sum, FUN = "/")
  eta_reformat<-matrix(0, n.cells, n.doc)
  
  plug<-function(target, id, data){
    target[id+1]<-data
    return(target)
  }
  
  eta_reformat<-mapply(plug, split(eta_reformat, row(eta_reformat)), 
             split(neigh_centers, row(neigh_centers)), 
             split(eta, row(eta)), SIMPLIFY = FALSE)
  eta_reformat<-do.call(rbind,eta_reformat)
  
  p_wk <- eta_reformat %*% theta_dk * beta_wk[C+1,]
  p_w <- apply(p_wk, 1, FUN = sum)
  log_pw <- log(p_w)
  res <- exp(-sum(log_pw)/n.cells)
  return(res)
}


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


## Let tissues is a list of data frames (X, Y, type)
gibbs_spatial_LDA_multiple<-function(tissue, ntopics, sigma, region_radius, kneigh,
                                     npoints_selected = 1, ini_LDA = 0, ninit = 1, 
                                     niter_init = 100, beta = .05, alpha = .01,
                                     trace = 0, compute_loglike = 0, seed = 123, thin = 20, burnin = 1000,
                                     niter = 200, save_image = 0, save_data = 0, display_progress = 1,
                                     output_path = "test/",fig.width = 30, fig.height = 30,
                                     do.parallel = 1, n.cores = 1){
  
  set.seed(seed)
  
  if(length(tissue) == 1) do.parallel<- 0 ## no need parallel computing if only one image
  
  ### combine multiple image into a single data frame
  itr_df<-do.call(rbind, tissue)
  itr_df$X<-as.numeric(itr_df$X)
  itr_df$Y<-as.numeric(itr_df$Y)
  itr_df$type<-as.factor(itr_df$type) 
  itr_df$image<-as.factor(itr_df$image)

  ## number of cells per image
  ncells<-table(itr_df$image)
  ## ncells<-unlist(lapply(tissue,nrow))
  print("number of cells per image:")
  print(ncells)
  
  ### for coords for each sample
  coords<-lapply(tissue,GetCoords)

  print("Initialization....")
  
  perplexity_min<-Inf
  Z_keep<-NULL
  D_keep<-NULL
  neigh_centers_keep<-NULL
  neigh_dists_keep<-NULL
  
  print("Numer of Initializations:")
  print(ninit)

  
  for (ini in 1:ninit){ ### should chose the best one among random sample points
    
  
  if(do.parallel){
    
    ### register cores for parallel computing 
    require(doParallel)
    registerDoParallel(n.cores)
  
    results<-foreach(i = 1:length(tissue),.packages=c('sp','RANN')) %dopar% {
      
      source("scripts/spatial_func.R")
      
      ## set seed for each core in parallel computing
      set.seed(seed+i)
      
      ## get coords
      coords_selected<-coords[[i]]
      coords_selected$X<-as.numeric(coords_selected$X)
      coords_selected$Y<-as.numeric(coords_selected$Y)
      
      ##spatial stratified sampling
      ## parameters need to be further selected 
      num_x_strata<-(max(coords_selected$X)-min(coords_selected$X))/(region_radius*2)
      num_y_strata<-(max(coords_selected$Y)-min(coords_selected$Y))/(region_radius*2)
      num_x_strata<-round(num_x_strata)
      num_y_strata<-round(num_y_strata)
      #print("Selected X, Y strata:\n")
      #print(num_x_strata)
      #print(num_y_strata)
      
      
      ### need a C function for efficiently spatial sampling (maybe rejective sampling)
      center_idx<-stratified_sampling_idx(coords_selected,npoints_selected,num_x_strata,num_y_strata)
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
    
    results<-foreach(i = 1:length(tissue),.packages=c('sp','RANN')) %do% {
      
      ## get coords
      coords_selected<-coords[[i]]
      coords_selected$X<-as.numeric(coords_selected$X)
      coords_selected$Y<-as.numeric(coords_selected$Y)
      
      ##spatial stratified sampling
      ## parameters need to be further selected 
      num_x_strata<-(max(coords_selected$X)-min(coords_selected$X))/(region_radius*2)
      num_y_strata<-(max(coords_selected$Y)-min(coords_selected$Y))/(region_radius*2)
      num_x_strata<-round(num_x_strata)
      num_y_strata<-round(num_y_strata)
      #print("Selected X, Y strata:\n")
      #print(num_x_strata)
      #print(num_y_strata)
      
      
      ### need a C function for efficiently spatial sampling (maybe rejective sampling)
      center_idx<-stratified_sampling_idx(coords_selected,npoints_selected,num_x_strata,num_y_strata)
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
  if(length(tissue) > 1){
    add<-c(0,unlist(lapply(1:(length(ncenters)-1),function(x) sum(ncenters[1:x]))))
    neigh_centers<-neigh_centers+rep(add,times = ncells)
  }
    
  ## Distance of cells to their closest centers 
  neigh_dists<-lapply(results, function(x) x$neigh_dists)
  neigh_dists<-do.call(rbind,neigh_dists)
  
  ############ randomly sample Z and D ---------------------------------------------------
  
  ## Initialize Z, D for spatial Topic model
  neigh_centers <- neigh_centers - 1L
  D<-as.integer(neigh_centers[,1])
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
  ## test Gibbs sampling codes for spatial LDA
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
    
    
    ### initialization with LDA
    if(ini_LDA){
      #system.time(gibbs.res <- gibbs_lda_c(
      #  docs,
      #  Ndk, Nwk, Nk, Z, D,
      #  doc_list,
      #  word_list,
      #  K, niter_init, beta,alpha))
      #perplexity<- perplexity(Ndk, Nwk, Nd, Nk, docs[,2], D, K, V, alpha = alpha, beta = beta)
      #print(perplexity)
      system.time(gibbs.res <-gibbs_sampler_c(docs, Ndk, Nwk, Nk, Nd, Z, D, neigh_dists, 
                                          neigh_centers, doc_list, word_list, 
                                          K, beta, alpha, sigma, thin, niter_init, 0,
                                          0, 0, 0))
      #perplexity<-perplexity_spaLDA(gibbs.res$Theta, gibbs.res$Beta, C, neigh_dists, sigma,
      #                                      neigh_centers, M, alpha, beta)
      #print(perplexity)
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
  }
  
  ### Keep the best results during initialization
  if(ini_LDA){
    print("Min perplexity during initialization:")
    print(perplexity_min)
  }
  Z<-Z_keep
  D<-D_keep
  neigh_centers<-neigh_centers_keep
  neigh_dists<-neigh_dists_keep

  print("number of centers selected:")
  print(ncenters) ### [TODO] change this later

  print("number of cells per region:")
  print(mean(table(D)))
  
  
  if(ninit > 1){ ## need update only when several initialization
      #K<-as.integer(ntopics)
      ## number of celltypes
      #V<-length(unique((C)))
      ## number of regions
      M<-length(table(D))  ## M might be changed since D changed
      ## number of words
      #N<-length(C)
    
    Ndk <- table_2d_fast(D, Z, M, K) ## number of cells per topic per region (image specific)
    Nwk <- table_2d_fast(C, Z, V, K) ## number of cells per topic per celltype
    Nk <- table_1d_fast(Z, K)  ## number of cells in each topic 
    Nd <- table_1d_fast(D, M) ## number of cells in each region (image specific)
    
    ## list of different docs/regions
    doc_list<-as.integer(1:M)-1L
    ## list of different words
    word_list<-as.integer(1:V)-1L
  }
  
  
  print("Finish initialization. Start Gibbs sampling....")
  
  
  ## Gibbs_sampling 
  ## [TODO] should allow multiple chains to be run in parallel.
  ## Simply add results from multiple chains
  ## [TODO] can use Rcppparallel
  system.time(gibbs.res <-gibbs_sampler_c(docs, Ndk, Nwk, Nk, Nd, Z, D, neigh_dists, 
                                          neigh_centers, doc_list, word_list, 
                                          K, beta, alpha, sigma, thin, burnin, niter,
                                          trace, display_progress, compute_loglike))
  
  
  ## calculate perplexity based on current document assignment
  #gibbs.res$perx<-perplexity(Ndk, Nwk, Nd, Nk, docs[,2], D, K, V, alpha = alpha, beta = beta)
  ##gibbs.res$perx<-perplexity_spaLDA(gibbs.res$Theta, gibbs.res$Beta, C, neigh_dists, sigma,
  ##                                            neigh_centers, M, alpha, beta)
  print("Output model perplexity..")
  print(gibbs.res$Perplexity)
  
  
  if(save_image){
    ## assign cell with the max prob
    prob<-as.matrix(gibbs.res$Z.trace)
    itr_df$Z_max<-as.factor(apply(prob,1,which.max))
    
    for(image in levels(itr_df$image)){
      
      ggplot(itr_df[itr_df$image == image, ],aes(X,Y,color = Z_max))+geom_point(size = 0.8)+
        scale_color_manual(limits = levels(itr_df$Z_max),values = hue_pal()(ntopics))
    
      ggsave(paste0(output_path,"test_K",K,"_sigma",sigma, "_radius",region_radius,"_center",
                  npoints_selected,"_kneigh",kneigh, "_iter",niter*thin,"_burn", burnin, 
                  "_initLDA",ini_LDA,"_",image,"_max_prob.png"),width = fig.width,height = fig.height,units = "in")
    }
  }
  
  if(save_data){
    save(docs, Ndk, Nwk, Nk, Nd,
         Z,D, neigh_dists, neigh_centers,
         doc_list, word_list,
         K, sigma, gibbs.res,
         file = paste0(output_path,"test_K",K,"_sigma",sigma, "_radius",region_radius,"_center",
                       npoints_selected,"_kneigh",kneigh,"_iter",niter*thin,"_burn",
                       burnin, "_initLDA",ini_LDA,"_max_prob.rdata"))
  }
  return(gibbs.res)
  
}

