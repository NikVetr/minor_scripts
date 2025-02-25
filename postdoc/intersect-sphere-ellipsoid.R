# -----------------------------------------------
# Function to solve for the intersection of an ellipse and a circle in 2D
# -----------------------------------------------
solve_ellipse_circle_2d <- function(A, B, r) {
  # A, B are the coefficients of the ellipse Ax^2 + By^2 = 1
  # r is the radius of the circle x^2 + y^2 = r^2
  # This function finds valid (x, y) pairs at the intersection
  
  f <- function(x) A * x^2 + B * (r^2 - x^2) - 1
  sol <- tryCatch(
    uniroot(f, lower = -r, upper = r, tol = 1e-9)$root,
    error = function(e) return(NULL)  # Catch numerical issues
  )
  
  if (is.null(sol)) return(NULL)  # Handle cases with no solution
  
  x_vals <- sol
  y_vals <- sqrt(r^2 - x_vals^2)  # Solve for y using the circle equation
  return(c(x_vals, y_vals))
}

# -----------------------------------------------
# Recursive function to sample from the intersection of hypersphere and hyperellipsoid
# -----------------------------------------------
recursive_sample <- function(L, hsrad = 1) {
  #
  # L: Vector of hyperellipsoid coefficients (e.g., c(2.5, 0.75, 0.5, 0.25) for 4D case)
  #
  # Output: A point in R^n that lies exactly on the intersection
  #
  
  n <- length(L)  # Dimension of the space
  if (n < 3) stop("Dimension must be at least 3 for recursion.")
  
  # Step 1: Find the largest coefficient in L and subtract using it
  which.med <- function(x) rank(x)[ceiling(length(x)/2)]
  med_index <- which.med(L)  # Find index of med coefficient
  L_med <- L[med_index]  # medimum coefficient
  L_remaining <- L[-med_index]  # All except the med one
  
  # Step 2: Compute derived ellipsoid coefficients (difference trick)
  M <- (L_med - L_remaining) / (L_med - hsrad)  # This ensures valid subtraction
  
  # Step 3: Compute minimum radius and valid range for r
  min_radius <- sqrt(1/max(M))  # Smallest axis-aligned distance from origin
  max_radius <- hsrad                 # Outer constraint from the unit sphere
  
  if (min_radius > max_radius) return(NULL)  # Handle invalid case
  
  # Step 4: Sample a radius r from a stretched Beta distribution
  r <- min_radius + (max_radius - min_radius) * rbeta(1, 2, 2)  # Adjust beta params as needed
  
  
  if (n == 3) {
    print("base case")
    # BASE CASE: We are now in 2D => Solve ellipse-circle intersection directly
    solution <- solve_ellipse_circle_2d(M[1], M[2], r)
    if (is.null(solution)) return(NULL)  # Catch cases with no valid solution
    return(c(solution, sqrt(1 - r^2)))  # Compute the final missing coordinate
  }
  
  # Step 5: Recursively sample from the (n-1)-dimensional intersection
  sub_point <- recursive_sample(M, r)  # Call the function with reduced dimension
  
  if (is.null(sub_point)) return(NULL)  # Handle failed recursive sampling
  
  # Step 6: Solve for the missing coordinate to ensure |x| = 1
  last_coord <- sqrt(1 - sum(sub_point^2))  # Solve for the last coordinate
  
  # Step 7: Return the full sampled point, restoring original order
  result <- numeric(n)
  result[med_index] <- last_coord  # Place the solved coordinate in its original position
  result[-med_index] <- sub_point  # Place the recursive result
  
  return(result)
}

# -----------------------------------------------
# Generate samples
# -----------------------------------------------
set.seed(42)
L <- c(2.5, 0.75, 0.5, 0.25)  # Example ellipsoid coefficients in 4D
points <- replicate(1000, recursive_sample(L), simplify = FALSE)
k <- do.call(rbind, points)  # Convert list of points to matrix

cat("Generated", nrow(points), "samples on the intersection.\n")

# Quick validation checks
valid_sphere <- sqrt(rowSums(points^2))  # Should be close to 1
valid_ellipsoid <- rowSums(sweep(points^2, 2, L, `*`))  # Should be close to 1

print(valid_sphere[1:10])
print(valid_ellipsoid[1:10])
