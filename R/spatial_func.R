## Some functions for spatial analysis

GetCoords <- function(tissue,axis = "2D") {
  if(axis == "3D"){
    coords <- as.data.frame(cbind(tissue$X,tissue$Y,tissue$Y2))
    colnames(coords) <- c("X", "Y", "Y2")
    return(coords)
  }else{
    coords <- as.data.frame(cbind(tissue$X,tissue$Y))
    colnames(coords) <- c("X", "Y")
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
  
  # List to store sampled points' indices
  sampled_point_idx <- list()
  
  points_in_cells<-t(st_within(points_sf, grid_sf))
  points_in_cells
  
  for(i in 1:length(grid_sf)){
    
    sampled_point_idx[[i]] <- ifelse(length(points_in_cells[[i]])==0, NA, 
                                   ifelse(length(points_in_cells[[i]]) < num_samples_per_stratum+1, points_in_cells[[i]], 
                                          sample(points_in_cells[[i]], num_samples_per_stratum)))
  }
  selected_points<- unlist(sampled_point_idx)
  selected_points <- selected_points[!is.na(selected_points)]
  
  return(selected_points)
}

