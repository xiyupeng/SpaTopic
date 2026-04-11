## Some functions for spatial analysis

GetCoords <- function(tissue,axis = "2D") {
  if(axis == "3D"){
    # Check for Z column first (preferred), then Y2 (legacy)
    if("Z" %in% colnames(tissue)){
      coords <- data.frame(
        X = as.numeric(tissue$X),
        Y = as.numeric(tissue$Y),
        Z = as.numeric(tissue$Z)
      )
    } else if("Y2" %in% colnames(tissue)){
      coords <- data.frame(
        X = as.numeric(tissue$X),
        Y = as.numeric(tissue$Y),
        Z = as.numeric(tissue$Y2)  # Rename Y2 to Z for consistency
      )
    } else {
      spatopic_message("ERROR", "For 3D analysis, tissue data must contain either 'Z' or 'Y2' column")
      return(NULL)
    }
    return(coords)
  }else{
    coords <- data.frame(
      X = as.numeric(tissue$X),
      Y = as.numeric(tissue$Y)
    )
    return(coords)
  }
}

#' Spatially stratified random sample points from an image.
#' 
#' Spatially stratified random sample points from an image by R package \code{sf}
#' 
#' @param points a data frame contains all points in a image with X, Y coordinates. 
#' 
#' @param cellsize a vector of length 2 contains the size of each grid square. 
#' Default c(600,600).
#' 
#' @param num_samples_per_stratum number of point selected from each grid square.
#' Default 1.
#' 
#' @return Return a vector contains index of sampled points.
#' 
#' @importFrom sf st_as_sf
#' @importFrom sf st_make_grid
#' @importFrom sf st_geometry
#' @importFrom sf st_cast
#' @importFrom sf st_within
#' @importFrom sf st_crs
#' 
#' @examples 
#' 
#' data("lung5")
#' pt_idx<-stratified_sampling_sf(lung5, cellsize = c(600,600))
#' 
#' @export
stratified_sampling_sf <- function(points, cellsize = c(600,600), num_samples_per_stratum = 1) {
  
  # Convert points to an sf object
  points_sf <- st_as_sf(points, coords = c("X", "Y"))
  
  # Create grid cells for stratification using cellsize
  grid_sf <- st_make_grid(points_sf, 
                          cellsize = cellsize, 
                          crs = st_crs(points_sf))
  grid_sf <- st_geometry(grid_sf)
  grid_sf <- st_cast(grid_sf, "POLYGON")
  
  points_in_cells<-t(st_within(points_sf, grid_sf))
  
  # List to store sampled points' indices
  sampled_point_idx <- lapply(points_in_cells, function(cell) {
    if (length(cell) == 0) {
      return(NA)
    } else if (length(cell) < num_samples_per_stratum + 1) {
      return(cell)
    } else {
      return(sample(cell, num_samples_per_stratum))
    }
  })
 
  
  selected_points<- unlist(sampled_point_idx)
  selected_points <- selected_points[!is.na(selected_points)]
  
  return(selected_points)
}

#' Spatially stratified random sample points from a 3D image using 2D slicing.
#' 
#' Spatially stratified random sample points from a 3D image by slicing along the Z-axis
#' and applying 2D stratified sampling using R package \code{sf} to each slice.
#' This approach leverages the efficiency of sf's 2D spatial indexing for 3D data.
#' 
#' @param points a data frame contains all points in a 3D image with X, Y, Z coordinates. 
#' 
#' @param cellsize a vector of length 2 contains the size of each grid square in X and Y dimensions. 
#' Default c(600,600).
#' 
#' @param z_cellsize numeric value specifying the thickness of each Z slice. 
#' Default 600.
#' 
#' @param num_samples_per_stratum number of points selected from each 3D grid cell (X,Y,Z).
#' Default 1.
#' 
#' @return Return a vector contains indices of sampled points from the original data frame.
#' 
#' @details
#' This function performs 3D stratified sampling by:
#' 1. Binning points into Z slices of thickness \code{z_cellsize}
#' 2. For each Z slice, applying 2D stratified sampling using \code{stratified_sampling_sf}
#' 3. Combining sampled indices from all slices
#' 
#' The approach is efficient for large 3D datasets as it leverages sf's optimized 
#' 2D spatial operations rather than implementing full 3D spatial indexing.
#' 
#' @seealso \code{\link{stratified_sampling_sf}} for 2D stratified sampling
#' 
#' @examples 
#' 
#' # Create example 3D data
#' set.seed(123)
#' tissue_3d <- data.frame(
#'   X = runif(1000, 0, 1000),
#'   Y = runif(1000, 0, 1000), 
#'   Z = runif(1000, 0, 100)
#' )
#' 
#' # Sample points with 3D stratification
#' pt_idx <- stratified_sampling_3D_via_2D(tissue_3d, 
#'                                          cellsize = c(200, 200), 
#'                                          z_cellsize = 200)
#' 
#' @export
stratified_sampling_3D_via_2D <- function(points, cellsize = c(600,600), z_cellsize = 600, num_samples_per_stratum = 1) {

  if (!all(c("X", "Y", "Z") %in% colnames(points))) {
    spatopic_message("ERROR", "Input data must have columns X, Y, and Z")
  }
  
  #points$X <- as.numeric(points$X)
  #points$Y <- as.numeric(points$Y)
  #points$Z <- as.numeric(points$Z)
  z_cellsize <- as.numeric(z_cellsize)
  
  if(is.na(z_cellsize) || z_cellsize <= 0) {
    spatopic_message("ERROR", "z_cellsize must be a positive numeric value")
    return(NULL)
  }
  
  points$orig_idx <- seq_len(nrow(points))  # Track original row indices
  
  # Bin Z into slices of thickness z_cellsize
  z_min <- min(points$Z)
  z_max <- max(points$Z)
  

  if(z_max - z_min < z_cellsize){
    sampled_idx <- stratified_sampling_sf(points = points, cellsize = cellsize, num_samples_per_stratum = num_samples_per_stratum)
 
  }else{

    z_bin <- floor((points$Z - z_min) / z_cellsize)
    points$z_bin <- z_bin
    z_bins <- unique(z_bin)
   
    sampled_idx <- NULL
    
    for (zb in z_bins) {
      slice <- points[points$z_bin == zb, ]
      if (nrow(slice) == 0) next
      idx_in_slice <- stratified_sampling_sf(slice, cellsize = cellsize, num_samples_per_stratum = num_samples_per_stratum)
      sampled_idx <- c(sampled_idx, slice$orig_idx[idx_in_slice])
    }
  }
  
  return(sampled_idx)
}

