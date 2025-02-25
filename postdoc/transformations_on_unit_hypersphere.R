# Install/load required packages
library(cmdstanr)
library(posterior)

#function
sample_sphere_normals <- function(k, n) {
  X <- matrix(rnorm(n * k), nrow = n, ncol = k)
  X <- X / sqrt(rowSums(X^2))
  X
}

# Write the Stan model as a string
stan_model_code <- "
data {
  // p is the ambient space dimension (the sphere is S^(p-1))
  int<lower=2> p;
}
parameters {
  // First p-2 angles are in [0, pi]
  vector<lower=0, upper=pi()>[p-2] theta1;
  // Last angle is in [0, 2pi)
  real<lower=0, upper=2*pi()> theta_last;
}
transformed parameters {
  // Concatenate the angles into a single vector theta of length p-1
  vector[p-1] theta;
  theta[1:(p-2)] = theta1;
  theta[p-1] = theta_last;
  
  // Transform spherical coordinates to Cartesian coordinates.
  // The formulas are:
  //   x1 = cos(theta[1])
  //   x2 = sin(theta[1])*cos(theta[2])
  //   ...
  //   x[p-1] = sin(theta[1])*...*sin(theta[p-2])*cos(theta[p-1])
  //   x[p]   = sin(theta[1])*...*sin(theta[p-2])*sin(theta[p-1])
  vector[p] x;
  x[1] = cos(theta[1]);
  real prod = sin(theta[1]);
  for (i in 2:(p-2)) {
    x[i] = prod * cos(theta[i]);
    prod = prod * sin(theta[i]);
  }
  x[p-1] = prod * cos(theta[p-1]);
  x[p]   = prod * sin(theta[p-1]);
}
model {
  // Add the Jacobian adjustments.
  // For a sphere in ℝ^p, for i = 1, ..., p-2, add (p-1-i)*log(sin(theta[i]))
  for (i in 1:(p-2))
    target += (p - 1 - i) * log(sin(theta[i]));
  // The last angle theta[p-1] is uniform on [0, 2pi), so no adjustment is needed.
}
generated quantities {
  // Output the Cartesian coordinates for inspection.
  vector[p] x_out = x;
}
"

# Save the Stan model to a file (this also returns the file path)
stan_file <- write_stan_file(stan_model_code)

# Compile the Stan model using cmdstanr
mod <- cmdstan_model(stan_file)

# Specify data: let's use p = 5 (so the sphere is S^4)
p <- 100
data_list <- list(p = p)

# Fit the model (adjust iterations/chains as needed)
fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Print a summary of the parameters and generated quantities
print(fit$summary(c("theta1", "theta_last", "x_out")))

# Extract the generated Cartesian coordinates from the draws
posterior_draws <- fit$draws("x_out", format = "df")

# Use bayesplot to compare the marginal distribution of x_out[1]

# Generate samples using both methods
n <- 1E3
samples_normals <- sample_sphere_normals(p, n)
samples_spherical <- as.data.frame(posterior_draws)[,1:p]

# Compare the marginal distributions.
# We overlay histograms for each coordinate.
par(mfrow = c(3, 5), mar = c(4, 4, 2, 1))
inds <- c(1, sample(1:p, min(13, p)), p)
for(i in inds) {
  hist(samples_normals[, i], breaks = -10:10/10, col = rgb(0, 0, 1, 0.5),
       main = paste("Coord", i, "\nNormals"),
       xlab = "", probability = TRUE, xlim = c(-1,1))
  hist(samples_spherical[, i], breaks = -10:10/10, col = rgb(1, 0, 0, 0.5),
       add = TRUE, probability = TRUE)
  legend("topright", legend = c("Normals", "Spherical"),
         fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)), bty = "n", cex = 0.8)
}

