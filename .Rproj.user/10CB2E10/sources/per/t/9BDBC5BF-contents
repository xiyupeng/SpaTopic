#library(sp)

GetCoords <- function(tissue) {
  coords <- as.data.frame(cbind(tissue$X,tissue$Y))
  colnames(coords) <- c("X", "Y")
  return(coords)
}

### need to rewrite it in C function
stratified_sampling_idx <- function(points, num_samples_per_stratum = 10, 
                                    num_x_strata = 5, num_y_strata =5) {
  # Convert points to a SpatialPointsDataFrame
  rownames(points)<-1:nrow(points)
  points_sp <- SpatialPoints(points)
  
  # Create vectors of stratum IDs based on X and Y coordinates
  x_strata <- cut(points[, "X"], breaks = num_x_strata, labels = FALSE)
  y_strata <- cut(points[, "Y"], breaks = num_y_strata, labels = FALSE)
  
  # Create an empty array to hold the indices of the sampled points
  sampled_point_idx <- list()
  
  # Sample points from each stratum
  for (i in 1:num_x_strata) {
    for (j in 1:num_y_strata) {
      stratum_points <- points_sp[x_strata == i & y_strata == j, ]
      num_points_in_stratum <- nrow(stratum_points@coords)
      if (num_points_in_stratum == 0) {
        next
      }
      if (num_points_in_stratum < num_samples_per_stratum) {
        warning(paste0("Stratum (", i, ",", j, ") has only ",
                       num_points_in_stratum, " points."))
        sampled_point_idx[[paste0(i,j)]] <- as.numeric(row.names(stratum_points))
      } else {
        sampled_point_idx[[paste0(i,j)]] <- as.numeric(sample(row.names(stratum_points),
                                               size = num_samples_per_stratum))
      }
    }
  }
  
  # Return the indices of the sampled points
  return(unlist(sampled_point_idx))
}

