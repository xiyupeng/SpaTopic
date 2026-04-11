#' SpaTopic: Spatial Topic Modeling for Multiplexed Images
#'
#' @description
#' SpaTopic is a package for spatial topic modeling of Multiplexed images.
#' It adapts an approach originally developed for image segmentation in computer vision, 
#' incorporating spatial information into the flexible design of regions (image partitions, 
#' analogous to documents in language modeling). We further refined the approach to address 
#' unique challenges in cellular images and provide an efficient C++ implementation of the 
#' algorithm in this R package.
#' 
#' Compared to other KNN-based methods (such as KNN-kmeans, the default neighborhood 
#' analysis in Seurat v5 R package), SpaTopic runs much faster on large-scale image datasets.
#'
#' The main functions in the 'SpaTopic' package
#' 
#' \itemize{
#' \item Prepare input   \code{\link{Seurat5obj_to_SpaTopic}}
#' \item Model Inference \code{\link{SpaTopic_inference}}
#' \item Print results \code{\link{print.SpaTopic}}
#' }
#'
#' @source <https://github.com/xiyupeng/SpaTopic>
#' @name SpaTopic-Package
#' 
#' @author Xiyu Peng \email{pansypeng124@gmail.com}
#'
#' @keywords package
#' 
#' @importFrom Rcpp sourceCpp
#' @importFrom stats var
#' @importFrom utils head
#' 
#' @useDynLib SpaTopic
"_PACKAGE"

