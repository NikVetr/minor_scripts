rotate_grid <- function(grid, theta) {
  theta_rad <- theta * (pi / 180)
  inv_theta_rad <- -theta_rad  # Inverse rotation
  
  # Dimensions of the original grid
  n <- nrow(grid)
  m <- ncol(grid)
  
  # Center of the original grid
  center_x <- (m + 1) / 2
  center_y <- (n + 1) / 2
  
  # Calculate extremes of the original grid's coordinates after rotation
  orig_coords <- expand.grid(x = c(1, m, 1, m), y = c(1, 1, n, n))
  rot_coords_x <- cos(theta_rad) * (orig_coords$x - center_x) - sin(theta_rad) * (orig_coords$y - center_y) + center_x
  rot_coords_y <- sin(theta_rad) * (orig_coords$x - center_x) + cos(theta_rad) * (orig_coords$y - center_y) + center_y
  x_extremes <- range(rot_coords_x)
  y_extremes <- range(rot_coords_y)
  
  # Dimensions of the new grid
  new_m <- ceiling(x_extremes[2] - x_extremes[1])
  new_n <- ceiling(y_extremes[2] - y_extremes[1])
  new_center_x <- (new_m + 1) / 2
  new_center_y <- (new_n + 1) / 2
  
  # New grid
  new_grid <- matrix(-2, nrow = new_n, ncol = new_m)
  
  # Apply inverse rotation to new grid coordinates
  new_coords <- expand.grid(x = 1:new_m, y = 1:new_n)
  x_inv_rotated <- round(cos(inv_theta_rad) * (new_coords$x - new_center_x) - sin(inv_theta_rad) * (new_coords$y - new_center_y) + center_x)
  y_inv_rotated <- round(sin(inv_theta_rad) * (new_coords$x - new_center_x) + cos(inv_theta_rad) * (new_coords$y - new_center_y) + center_y)
  
  # Valid indices within the original grid
  valid <- x_inv_rotated >= 1 & x_inv_rotated <= m & y_inv_rotated >= 1 & y_inv_rotated <= n
  
  # Map valid coordinates from original to new grid
  new_grid[cbind(new_coords$y[valid], new_coords$x[valid])] <- grid[cbind(y_inv_rotated[valid], x_inv_rotated[valid])]
  
  return(new_grid)
}

unrotate_grid <- function(orig_grid, rotated_grid, theta, new_values) {
  # Reverse rotation
  unrotated_grid <- rotate_grid(rotated_grid, -theta)
  
  # Calculate dimensions and offsets
  orig_n <- nrow(orig_grid)
  orig_m <- ncol(orig_grid)
  unrot_n <- nrow(unrotated_grid)
  unrot_m <- ncol(unrotated_grid)
  offset_row <- floor((unrot_n - orig_n) / 2) + 1
  offset_col <- floor((unrot_m - orig_m) / 2) + 1
  
  # Trim to match dimensions of original grid
  trimmed_unrotated_grid <- unrotated_grid[offset_row:(offset_row + orig_n - 1), offset_col:(offset_col + orig_m - 1)]
  
  # Update original grid with new values where original grid has 0s
  update_mask <- orig_grid == 0 & trimmed_unrotated_grid %in% new_values
  orig_grid[update_mask] <- trimmed_unrotated_grid[update_mask]
  
  return(orig_grid)
}

# Example usage
orig_grid <- matrix(0, 10, 10)
orig_grid[4:6, ] <- 1
theta <- 45
rotated_grid <- rotate_grid(orig_grid, theta)
# Assume modifications are made to rotated_grid here

rotated_grid[rotated_grid == 0] <- 3
new_values <- c(2, 3, 4)  # Example new values

# Unrotate the grid
unrotated_modified_grid <- unrotate_grid(orig_grid, rotated_grid, theta, new_values)
print(unrotated_modified_grid)



#### expand grid ####

neighbors <- function(grid, pts, val = 0, distance = 1){
  dims <- dim(grid)
  offset_grid <- expand.grid(y = -distance:distance, 
                             x = -distance:distance)
  offset_grid <- offset_grid[!apply(offset_grid, 1, function(x) all(x==0)),]
  neighboring_pts <- pts[rep(1:nrow(pts), each = nrow(offset_grid)), ] + offset_grid[rep(1:nrow(offset_grid), times = nrow(pts)), ]
  neighboring_pts <- neighboring_pts[!duplicated(neighboring_pts),]
  neighboring_pts <- neighboring_pts[neighboring_pts[,1] >= 1 & neighboring_pts[,1] <= dims[1] &
                                       neighboring_pts[,2] >= 1 & neighboring_pts[,2] <= dims[2],]
  return(as.matrix(neighboring_pts[grid[as.matrix(neighboring_pts)] == val,]))
}


has_neighbors <- function(grid, pts, val = 0, distance = 1){
  
  #make useful variables
  dims <- dim(grid)
  offset_grid <- expand.grid(y = -distance:distance, 
                             x = -distance:distance)
  offset_grid <- offset_grid[!apply(offset_grid, 1, function(x) all(x==0)),]
  
  #find neighboring points
  pt_inds <- rep(1:nrow(pts), each = nrow(offset_grid))
  neighboring_pts <- pts[pt_inds, ] + offset_grid[rep(1:nrow(offset_grid), times = nrow(pts)),]
  
  #subset to only valid (inside of grid) points
  in_grid <- neighboring_pts[,1] >= 1 & neighboring_pts[,1] <= dims[1] &
    neighboring_pts[,2] >= 1 & neighboring_pts[,2] <= dims[2]
  pt_inds <- pt_inds[in_grid]
  neighboring_pts <- neighboring_pts[in_grid,]
  
  #find desired points
  target_pt_inds <- unique(pt_inds[grid[as.matrix(neighboring_pts)] == val])
  
  #return pt values
  return(as.matrix(pts[target_pt_inds,]))
}

grid <- matrix(0, 10, 10)
grid[4:6, 2:3] <- 1
grid[7:9, 2:5] <- 2
grid[2:5, 4:5] <- 3
plot(which(grid == 1, arr.ind = T), xlim = c(-5, 15), ylim = c(-5, 15), pch = 19, col = 1)
points(which(grid == 2, arr.ind = T), pch = 19, col = 2)
points(which(grid == 3, arr.ind = T), pch = 19, col = 3)
points(which(grid == 0, arr.ind = T), pch = 19, col = "grey90")
points(neighbors(grid, pts = which(grid == 1, arr.ind = T), val = 0, distance = 1), pch = 19, col = 4)

expand_blocks <- function(grid, pt_sets = NULL, prop_expansion_rate = F, max_expansion_rate = 3) {
  
  if(is.null(pt_sets)){
    grid_vals <- sort(unique(c(0, grid)))
    pt_sets <- lapply(grid_vals, function(i){
      has_neighbors(grid, which(grid == i, T), 0)
    })
    names(pt_sets) <- grid_vals
  }
  
  if(prop_expansion_rate){
    n_border_pts <- sapply(pt_sets[as.character(grid_vals[grid_vals != 0])], nrow)
    expansion_distances <- ceiling(n_border_pts / max(n_border_pts) * max_expansion_rate)
  } else {
    expansion_distances <- setNames(rep(1, length(grid_vals) - 1), grid_vals[grid_vals != 0])
  }
  
  exp_grid <- grid
  while(nrow(pt_sets[["0"]]) != 0){
    expansion_order <- sample(setdiff(grid_vals, 0))
    for(i in expansion_order){
      
      if(nrow(pt_sets[[as.character(i)]]) == 0){
        next()
      }
      
      #find zero-valued neighbors
      neighboring_pts <- neighbors(grid = exp_grid,
                                   pts = pt_sets[[as.character(i)]],
                                   val = 0,
                                   distance = expansion_distances[as.character(i)])
      
      if(nrow(neighboring_pts) == 0){
        pt_sets[[as.character(i)]] <- neighboring_pts
        next()
      } else {
        #swap in to pt_sets and mark in grid
        exp_grid[neighboring_pts] <- i
        pt_sets[[as.character(i)]] <- has_neighbors(exp_grid, neighboring_pts, 0)
        
        #remove from 0-valued pt_set
        zero_vals_comb <- rbind(neighboring_pts, pt_sets[["0"]])
        pt_sets[["0"]] <- pt_sets[["0"]][!duplicated(zero_vals_comb)[-(1:nrow(neighboring_pts))],]
      }
      
    }
  }
  
  return(exp_grid)
}

# Example usage
grid <- matrix(0, 10, 10)
grid[4:6, 2:3] <- 1
grid[7:9, 2:5] <- 2
grid[2:5, 4:5] <- 3
plot(which(grid == 1, arr.ind = T), xlim = c(-5, 15), ylim = c(-5, 15), pch = 19, col = 1)
points(which(grid == 2, arr.ind = T), pch = 19, col = 2)
points(which(grid == 3, arr.ind = T), pch = 19, col = 3)
points(which(grid == 0, arr.ind = T), pch = 19, col = "grey90")

expanded_grid <- expand_blocks(grid, prop_expansion_rate = T)
plot(which(expanded_grid == 1, arr.ind = T), xlim = c(-5, 15), ylim = c(-5, 15), pch = 19, col = 1)
points(which(expanded_grid == 2, arr.ind = T), pch = 19, col = 2)
points(which(expanded_grid == 3, arr.ind = T), pch = 19, col = 3)
points(which(expanded_grid == 0, arr.ind = T), pch = 19, col = "grey90")


grid <- matrix(0, 100, 100)
grid[45:55, 45:55] <- 1  # Example block
grid[80:90, 80:90] <- 2  # Another block
grid[11:60, 10:12] <- 3  # Another block

plot(which(grid == 1, arr.ind = T), xlim = c(-5, 105), ylim = c(-5, 105), pch = 19, col = 1)
points(which(grid == 2, arr.ind = T), pch = 19, col = 2)
points(which(grid == 3, arr.ind = T), pch = 19, col = 3)
points(which(grid == 0, arr.ind = T), pch = 19, col = "grey90")

expanded_grid <- expand_blocks(grid, prop_expansion_rate = F)
plot(which(expanded_grid == 1, arr.ind = T), xlim = c(-5, 105), ylim = c(-5, 105), pch = 19, col = 1)
points(which(expanded_grid == 2, arr.ind = T), pch = 19, col = 2)
points(which(expanded_grid == 3, arr.ind = T), pch = 19, col = 3)
points(which(expanded_grid == 0, arr.ind = T), pch = 19, col = "grey90")
