# Load required libraries
library(kernlab)
library(scatterplot3d)
library(stats)
library(ape)
library(phangorn)
library(destiny)
library(vegan)
library(torch)

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
num_points <- 200
k <- 2  # User-specified dimension of the sphere
r1 <- 1 #radius of the first hypersphere
r2 <- 3 #squared radius of the second hypersphere

# Generate points on the surface of the first k-sphere centered at the origin
points1 <- sample_ksphere(num_points, k) * r1

# Generate points on the surface of the second k-sphere centered at (1, 0, 0, ..., 0) (shift along the first axis)
points2 <- sample_ksphere(num_points, k) * r2

#or from some other surface
get_shape <- function(p, n = 1000) {
  # Define theta for parametric representation
  # theta <- seq(0, 2 * pi, length.out = n)
  theta <- runif(n, 0, 2*pi)
  
  # Compute x and y coordinates
  x <- cos(theta)^p
  y <- sin(theta)^p
  
  return(data.frame(x=x, y=y))
}
p <- 3
points2 <- as.matrix(get_shape(p, n = num_points))
points2 <- points2 * r2

target_overlap <- 0.0
disp <- 0.0
# disp <- find_displacement(points1, points2, p = target_overlap, max_distance = 20)
points2[, 1] <- points2[, 1] + disp  # Shift center by 1 along the first axis

# Combine points into a single dataset
points <- rbind(points1, points2)
labels <- factor(c(rep(1, num_points), rep(2, num_points)))  # Labels for visualization

#convert to spherical?
cartesian_to_spherical <- function(cart) {
  cart <- as.matrix(cart)
  n <- nrow(cart)
  p <- ncol(cart)
  
  r <- sqrt(rowSums(cart^2))
  theta <- matrix(NA, nrow = n, ncol = p - 1)
  
  for (k in 1:(p - 1)) {
    denom <- sqrt(rowSums(cart[, k:p, drop = FALSE]^2))
    
    # Handle boundary case
    theta_k <- rep(0, n)
    nonzero <- denom > 0
    
    if (k == p - 1) {
      theta_k[nonzero] <- atan2(cart[nonzero, p], cart[nonzero, p - 1])
    } else {
      theta_k[nonzero] <- acos(cart[nonzero, k] / denom[nonzero])
    }
    
    theta[, k] <- theta_k
  }
  
  colnames(theta) <- paste0("theta", 1:(p - 1))
  result <- cbind(r = r, theta)
  return(result)
}

# points <- cartesian_to_spherical(points)

# Step 2: Apply kPCA with different kernels
# RBF Kernel
kpca_rbf <- kpca(~., data = as.data.frame(points), kernel = "rbfdot", 
                 kpar = list(sigma = 0.5), features = 4)
points_rbf_transformed <- rotated(kpca_rbf)
eig(kpca_rbf) / sum(eig(kpca_rbf))

# other kernel
kern2use <- c("rbfdot", #1
              "polydot", #2
              "vanilladot", #3
              "tanhdot", #4
              "laplacedot", #5
              "besseldot",#6
              "anovadot",#7
              "splinedot")[c(8)]

points_other_transformed <- points
for(i in 1:length(kern2use)){
  kpca_other <- kpca(~., data = as.data.frame(points_other_transformed), kernel = kern2use[i], features = 2, kpar = list())
  points_other_transformed <- rotated(kpca_other)
}

# Linear Kernel
kpca_linear <- prcomp(points)
points_linear_transformed <- kpca_linear$x[,1:2]

#laplacian diffusion map
diff_map <- DiffusionMap(points)
diff_coords <- eigenvectors(diff_map)

#try isomap out
# dist_matrix <- dist(points)
# isomap_result <- isomap(dist_matrix, k = k, ndim = 2)
# isomap_coords <- scores(isomap_result)
isomap_coords <- NA

#try umap and tsne
umap_result <- uwot::umap(points)
tsne_result <- Rtsne::Rtsne(points)$Y

#try a small VAE
points_tensor <- torch_tensor(as.matrix(points))

VAE <- nn_module(
  "VAE",
  
  initialize = function(input_dim, hidden_dims = c(64, 32), latent_dim = 2) {
    # Store input dimensions
    self$input_dim <- input_dim
    self$latent_dim <- latent_dim
    self$hidden_dims <- hidden_dims
    
    # -------------------
    # Define Encoder
    # -------------------
    self$enc_fc1 <- nn_linear(input_dim, hidden_dims[1])
    self$enc_fc2 <- nn_linear(hidden_dims[1], hidden_dims[2])
    self$enc_mu <- nn_linear(hidden_dims[2], latent_dim)        # Mean output
    self$enc_logvar <- nn_linear(hidden_dims[2], latent_dim)    # Log variance output
    
    # -------------------
    # Define Decoder
    # -------------------
    self$dec_fc1 <- nn_linear(latent_dim, hidden_dims[2])
    self$dec_fc2 <- nn_linear(hidden_dims[2], hidden_dims[1])
    self$dec_out <- nn_linear(hidden_dims[1], input_dim)        # Reconstructed output
  },
  
  encode = function(x) {
    # Forward pass through encoder
    h1 <- torch_relu(self$enc_fc1(x))
    h2 <- torch_relu(self$enc_fc2(h1))
    mu <- self$enc_mu(h2)
    logvar <- self$enc_logvar(h2)
    list(mu, logvar)
  },
  
  reparameterize = function(mu, logvar) {
    std <- torch_exp(0.5 * logvar)
    eps <- torch_randn_like(std)
    mu + eps * std  # Reparameterized latent vector
  },
  
  decode = function(z) {
    # Forward pass through decoder
    h1 <- torch_relu(self$dec_fc1(z))
    h2 <- torch_relu(self$dec_fc2(h1))
    torch_sigmoid(self$dec_out(h2))  # Sigmoid activation for output
  },
  
  forward = function(x) {
    # Full VAE pipeline: Encoding -> Reparameterization -> Decoding
    encoded <- self$encode(x)
    z <- self$reparameterize(encoded[[1]], encoded[[2]])
    decoded <- self$decode(z)
    list(decoded, encoded[[1]], encoded[[2]])  # Return decoded output, mu, logvar
  }
)

#loss function
vae_loss <- function(recon_x, x, mu, logvar) {
  recon_loss <- nnf_mse_loss(recon_x, x, reduction = "sum")
  kl_loss <- -0.5 * torch_sum(1 + logvar - mu^2 - torch_exp(logvar))
  recon_loss + kl_loss
}

#train VAE

# Define model parameters
hidden_dims <- c(8, 4)
latent_dim <- k

# Instantiate VAE model
vae <- VAE(k, hidden_dims, latent_dim)
optimizer <- optim_adam(vae$parameters, lr = 1e-3)

# Training loop
epochs <- 1
for (epoch in seq_len(epochs)) {
  optimizer$zero_grad()
  
  output <- vae(points_tensor)
  
  loss <- vae_loss(output[[1]], points_tensor, output[[2]], output[[3]])
  loss$backward()
  optimizer$step()
  
  if (epoch %% 10 == 0) {
    cat(sprintf("Epoch %d, Loss: %.2f\n", epoch, loss$item()))
  }
}

VAE_repr <- as.matrix(vae$encode(points_tensor)[[1]])

# Step 3: Hierarchical clustering
# dist_pts <- dist(points)
# hclust_result <- hclust(dist_pts)
# clusters <- cutree(hclust_result, k = 2)
# 
# nj_result <- nj(dist_pts)
# upgma_result <- upgma(dist_pts)
# clusters_upgma <- cutree(as.hclust(upgma_result), k = 2)

#### plot results ####
par(mfrow = c(3,4), mar = c(2,2,1,1))
if (k == 2) {
  # Plot original points in 2D if k = 2
  # par(mfrow = c(2,2), mar = c(4,4,2,2))
  plot(points[,1], points[,2], col = as.numeric(labels), 
       pch = 19, main = paste0("Orig. Pts (n = ", num_points, ")"), 
       xlab = "", ylab = "")
}

if (k == 3) {
  # Plot original points in 3D if k = 3
  # par(mfrow = c(1,1), mar = c(4,4,2,2))
  scatterplot3d(points[,1], points[,2], points[,3], color = as.numeric(labels), pch = 19, main = "Original Points in 3D", xlab = "X", ylab = "Y", zlab = "Z")
}

# Plot transformed points by kPCA (RBF)
plot(points_rbf_transformed[,1:2], col = as.numeric(labels), pch = 19, main = "kPCA (RBF)", xlab = "", ylab = "")

# Plot transformed points by kPCA (RBF) later PCs
# plot(points_rbf_transformed[,3:4], col = as.numeric(labels), pch = 19, main = "kPCA (RBF)", xlab = "PC3", ylab = "PC4")
# legend("topright", legend = c("hsph 1", "hsph 2"), col = c(1, 2), pch = 19)

# Plot transformed points by kPCA (othernomial)
plot(points_other_transformed, col = as.numeric(labels), pch = 19, main = paste0("kPCA (", kern2use, ")"), xlab = "", ylab = "")

# Plot transformed points by kPCA (Linear)
plot(points_linear_transformed, col = as.numeric(labels), pch = 19, main = "Linear PCA", xlab = "", ylab = "")

# Plot transformed points by diffusion map
plot(diff_coords, col = as.numeric(labels), pch = 19, main = "Diffusion Mapping", xlab = "", ylab = "")

# Plot transformed points by isomap
if(!is.na(isomap_coords)){
  plot(isomap_coords, col = as.numeric(labels), pch = 19, main = "Isomapping", xlab = "", ylab = "")
}

# Plot transformed points by VAE
for(i in 1:floor(ncol(VAE_repr) / 2)){
  plot(VAE_repr[,i*2-1:0], col = adjustcolor(as.numeric(labels), 1), pch = 19, main = "VAE bn", 
       xlab = paste0("Axis ", i*2-1), ylab = paste0("Axis ", i*2))
}


# Plot transformed points by umap
plot(umap_result, col = as.numeric(labels), pch = 19, main = "UMAP", xlab = "", ylab = "")

#plot tSNE results
plot(tsne_result, col = as.numeric(labels), pch = 19, main = "tSNE", xlab = "", ylab = "")


plot.new()
legend("topright", legend = c("hypersphere 1", "hypersphere 2"), col = c(1, 2), pch = 19)

#clusters obtained from hierarchical clustering
# mean(clusters == labels)
# mean(clusters_upgma == labels)

