## prepare input for SpaTopic

#' Convert a Seurat v5 object as the input of SpaTopic
#' 
#' Prepare SpaTopic input from one Seurat v5 object
#' 
#' @param object Seurat v5 object
#' 
#' @param group.by \code{character}. The name of the column that contains celltype information in the Seurat object.
#' 
#' @param image \code{character}. The name of the image. Default is "image1".
#' 
#' @return Return a data frame as the input of SpaTopic
#' 
#' @examples 
#' 
#' ## nano.obj is a Seurat v5 object
#' #dataset<-Seurat5obj_to_SpaTopic(object = nano.obj, 
#' #                 group.by = "predicted.annotation.l1",image = "image1")
#' ## Expect output
#' data("lung5")
#' 
#' @seealso \code{\link{lung5}}
#' 
#' @export
Seurat5obj_to_SpaTopic<-function(object,group.by,image = "image1"){
  
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    warning("The SeuratObject package (>= 5.0) must be installed to use this functionality")
    return(NULL)
  } else {
    
    if(is.list(object)){
      warning("The input is a list, not a seurat object. Just grab the first item in the list")
      object<-object[[1]]
    }
    
    tryCatch(
    {coords <- SeuratObject::GetTissueCoordinates(object, which = "centroids")},
    error = function(e){
      warning("Please check if the input object is a Seurat v5 object")
      return(NULL)
    },
    warning = function(w){
      print(w)
    }
    )
    
    
    tryCatch(
      {celltype<-object[[group.by]]},
      error = function(e){
        print(e)
        return(NULL)
      },
      warning = function(w){
        print(w)
      }
    )
    
    cells <- coords$cell
    rownames(coords) <- cells
    coords <- as.matrix(coords[, c("x", "y")])
    data_select<-as.data.frame(cbind(image, coords[,1], coords[,2], celltype))
    colnames(data_select)<-c("image","X","Y","type")
    data_select$image<-as.factor(data_select$image)
    data_select$type<-as.factor(data_select$type)
    return(data_select)
  }
}
