# Function to compute the analytical distance d
# Function to compute the overlap fraction given a displacement d
overlap_fraction <- function(d, x, y) {
  # Shift the second hypersphere by distance d along the first axis
  y_shifted <- y
  y_shifted[, 1] <- y_shifted[, 1] + d
  
  # Calculate distances from the origin for shifted points
  distances <- sqrt(rowSums(y_shifted^2))
  
  # Calculate the fraction of points within the original hypersphere
  overlap_frac <- mean(distances <= 1)
  return(overlap_frac)
}

# Function to find the displacement distance for a specified overlap fraction p
find_displacement <- function(x, y, p, max_distance = 3) {
  # Objective function: difference between current and desired overlap fraction
  f <- function(d) {
    frac <- overlap_fraction(d, x, y)
    return(frac - p)
  }
  
  # Ensure that the function crosses zero in the interval [lower, upper]
  lower <- 0
  upper <- 0.5  # Initial guess for upper bound
  while (f(upper) > 0 && upper < max_distance) {
    upper <- upper + 0.5
  }
  
  if (upper >= max_distance) {
    stop("Cannot find an upper bound where the overlap fraction is less than the desired p.")
  }
  
  # Use uniroot to find the displacement where f(d) = 0
  result <- uniroot(f, lower = lower, upper = upper)
  return(result$root)
}


# Function to generate random points within a unit hypersphere
generate_points_in_hypersphere <- function(n, num_points) {
  # Generate random points from a normal distribution
  points <- matrix(rnorm(n * num_points), nrow = num_points, ncol = n)
  # Normalize to lie within the unit hypersphere
  radii <- runif(num_points)^(1 / n)
  norms <- sqrt(rowSums(points^2))
  points <- (points / norms) * radii
  return(points)
}

# Simulation function
simulate_overlap <- function(n, num_points = 1e3) {
  
  # Generate points in the first hypersphere centered at the origin
  points1 <- generate_points_in_hypersphere(n, num_points)
  
  # Generate points in the second hypersphere centered at (d, 0, ..., 0)
  points2 <- generate_points_in_hypersphere(n, num_points)
  
  d <- find_displacement(points1, points2, 0.5)
  
  points2[, 1] <- points2[, 1] + d  # Shift along the first axis
  
  # Count the number of points from points2 that are within the original hypersphere
  distances <- sqrt(rowSums(points2^2))
  num_inside <- sum(distances <= 1)
  
  # Calculate the proportion of overlap
  overlap_fraction <- num_inside / num_points
  
  return(list(
    Dimension = n,
    Distance = d,
    OverlapFraction = overlap_fraction
  ))
}

# Run the simulation for dimensions 1 to 10
set.seed(42)  # For reproducibility
dimensions <- 1:10
num_points <- 1e3  # Number of points to sample
results <- lapply(dimensions, function(n) simulate_overlap(n, num_points))

# Convert results to a data frame
results_df <- do.call(rbind, lapply(results, as.data.frame))
print(results_df)
