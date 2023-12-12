#' Example input data for SpaTopic
#'
#' A lung cancer tumor tissue image, as the input of SpaTopic
#'
#' @format ## `lung5`
#' A data frame with 100149 rows and 4 columns: 
#' \describe{ 
#'   \item{image}{Image ID}
#'   \item{X}{X coordinate of the cell}
#'   \item{Y}{Y coordinate of the cell}
#'   \item{type}{Cell type}
#' } 
#' @source <https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/>
#' @seealso \code{\link{SpaTopic_inference}},\code{\link{Seurat5obj_to_SpaTopic}}
"lung5"
