#library(sp)
#library(sf)

GetCoords <- function(tissue,axis = "2D") {
  if(axis == "3D"){
    coords <- as.data.frame(cbind(tissue$X,tissue$Y,tissue$Y2))
    colnames(coords) <- c("X", "Y","Y2")
    return(coords)
  }else{
    coords <- as.data.frame(cbind(tissue$X,tissue$Y))
    colnames(coords) <- c("X", "Y")
    return(coords)
  }
}


stratified_sampling_idx_sf <- function(points, cellsize = c(600,600),num_samples_per_stratum = 1) {
  
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


