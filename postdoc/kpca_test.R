# Load required libraries
library(kernlab)
library(scatterplot3d)
library(stats)
library(ape)
library(phangorn)

#functions

# Function to compute the overlap fraction given a displacement d
overlap_fraction <- function(d, x, y) {
  y_shifted <- y
  y_shifted[, 1] <- y_shifted[, 1] + d
  distances <- sqrt(rowSums(y_shifted^2))
  overlap_frac <- mean(distances <= 1)
  return(overlap_frac)
}

# Function to find the displacement distance for a specified overlap fraction p
find_displacement <- function(x, y, p, max_distance = 3) {

  f <- function(d) {
    frac <- overlap_fraction(d, x, y)
    return(frac - p)
  }
  
  # Ensure that the function crosses zero in the interval [lower, upper]
  lower <- 0
  upper <- max_distance / 2  # Initial guess for upper bound
  while (f(upper) > 0 && upper < max_distance) {
    upper <- upper + 0.5
  }
  
  # Use uniroot to find the displacement where f(d) = 0
  result <- uniroot(f, lower = lower, upper = upper)
  return(result$root)
}


# Step 1: Function to generate points on the surface of a k-sphere
sample_ksphere <- function(num_points, k) {
  # Generate k standard normal random variables for each point
  points <- matrix(rnorm(num_points * k), ncol = k)
  # Rescale each point to unit length
  points <- t(apply(points, 1, function(row) row / sqrt(sum(row^2))))
  return(points)
}

# Number of points to sample per hypersphere
num_points <- 300
k <- 10  # User-specified dimension of the sphere

# Generate points on the surface of the first k-sphere centered at the origin
points1 <- sample_ksphere(num_points, k)

# Generate points on the surface of the second k-sphere centered at (1, 0, 0, ..., 0) (shift along the first axis)
points2 <- sample_ksphere(num_points, k) * 1.5

target_overlap <- 0.4
disp <- 0
# disp <- find_displacement(points1, points2, p = target_overlap, max_distance = 20)
points2[, 1] <- points2[, 1] + disp  # Shift center by 1 along the first axis

# Combine points into a single dataset
points <- rbind(points1, points2)
labels <- factor(c(rep(1, num_points), rep(2, num_points)))  # Labels for visualization

# Step 2: Apply Kernel PCA with different kernels
# RBF Kernel
kpca_rbf <- kpca(~., data = as.data.frame(points), kernel = "rbfdot", kpar = list(sigma = 0.5), features = 2)
points_rbf_transformed <- rotated(kpca_rbf)

# other kernel
kern2use <- "splinedot"
kpca_other <- kpca(~., data = as.data.frame(points), kernel = kern2use, features = 2, kpar = list())
points_other_transformed <- rotated(kpca_other)

# Linear Kernel
kpca_linear <- prcomp(points)
points_linear_transformed <- kpca_linear$x[,1:2]

# Step 3: Hierarchical clustering
dist_pts <- dist(points)
hclust_result <- hclust(dist_pts)
clusters <- cutree(hclust_result, k = 2)

nj_result <- nj(dist_pts)
upgma_result <- upgma(dist_pts)
clusters_upgma <- cutree(as.hclust(upgma_result), k = 2)

# Step 4: Visualize the original points and the Kernel PCA results
par(mfrow = c(2,2), mar = c(4,4,2,2))
if (k == 2) {
  # Plot original points in 2D if k = 2
  par(mfrow = c(2,2), mar = c(4,4,2,2))
  plot(points[,1], points[,2], col = as.numeric(labels), 
       pch = 19, main = "Original Points in 2D", 
       xlab = "X", ylab = "Y")
}

if (k == 3) {
  # Plot original points in 3D if k = 3
  par(mfrow = c(1,1), mar = c(4,4,2,2))
  scatterplot3d(points[,1], points[,2], points[,3], color = as.numeric(labels), pch = 19, main = "Original Points in 3D", xlab = "X", ylab = "Y", zlab = "Z")
}

# Plot transformed points by Kernel PCA (RBF)
plot(points_rbf_transformed, col = as.numeric(labels), pch = 19, main = "Points after Kernel PCA (RBF)", xlab = "PC1", ylab = "PC2")
legend("topright", legend = c("Hypersphere 1", "Hypersphere 2"), col = c(1, 2), pch = 19)

# Plot transformed points by Kernel PCA (othernomial)
plot(points_other_transformed, col = as.numeric(labels), pch = 19, main = paste0("Points after Kernel PCA (", kern2use, ")"), xlab = "PC1", ylab = "PC2")
```legend("topright", legend = c("Hypersphere 1", "Hypersphere 2"), col = c(1, 2), pch = 19)

# Plot transformed points by Kernel PCA (Linear)
plot(points_linear_transformed, col = as.numeric(labels), pch = 19, main = "Points after Kernel PCA (Linear)", xlab = "PC1", ylab = "PC2")
legend("topright", legend = c("Hypersphere 1", "Hypersphere 2"), col = c(1, 2), pch = 19)

#clusters obtained from hierarchical clustering
mean(clusters == labels)
mean(clusters_upgma == labels)
