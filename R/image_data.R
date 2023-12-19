#' Example input data for SpaTopic
#'
#' multiplexed image data on tumor tissue sample from non small cell lung cancer patient
#'
#' @format ## `lung5`
#' A data frame with 100149 rows and 4 columns: 
#' \describe{ 
#'   \item{image} Image ID
#'   \item{X} X coordinate of the cell
#'   \item{Y} Y coordinate of the cell
#'   \item{type} cell type
#' } 
#' @source <https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/>
#' @seealso \code{\link{SpaTopic_inference}},\code{\link{Seurat5obj_to_SpaTopic}}
"lung5"
