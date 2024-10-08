#' 'SpaTopic' R package 
#'
#' The 'SpaTopic' R package is centered around the 'SpaTopic' algorithm to infer the spatial
#' tissue architectures from multiplexed images. 
#' 
#' The package implements a Collapsed 
#' Gibbs sampling algorithm to infer topics, corresponding to distinct tissue 
#' microenvironments across multiple tissue images. 
#' Without obtaining the cell neighborhood info of every single cell, 
#' 'SpaTopic' runs much faster than other KNN-based methods on large-scale images.
#'
#' The main functions in the 'SpaTopic' package
#' 
#' \itemize{
#' \item Prepare input   \code{\link{Seurat5obj_to_SpaTopic}}
#' \item Model Inference \code{\link{SpaTopic_inference}}
#' }
#'
#'
#' @name SpaTopic-Package
#' 
#' @author Xiyu Peng \email{pansypeng124@gmail.com}
#'
#' @keywords package
#' 
#' @importFrom Rcpp sourceCpp
#' 
#' @useDynLib SpaTopic
"_PACKAGE"

