library(dplyr)
library(purrr)
library(parallel)

# Define the function to compute the inverse CDF
compute_inverse_cdf <- function(theta_q, w1, lambda, sigma) {
  # Constraints
  w2 <- w3 <- (1 - w1) / 2
  
  # Compute the inverse CDF based on theta_q and mixture parameters
  if (theta_q < w2) {
    # Left truncated exponential component
    u <- theta_q / w2
    Z_E <- 1 - exp(-lambda)
    x <- - (1 / lambda) * log(1 - u * Z_E)
  } else if (theta_q < (w2 + w1)) {
    # Middle truncated logistic component
    u <- (theta_q - w2) / w1
    Logistic_L <- plogis((0 - 0.5) / sigma)
    Logistic_U <- plogis((1 - 0.5) / sigma)
    inv_Logistic_arg <- u * (Logistic_U - Logistic_L) + Logistic_L
    inv_Logistic_arg <- pmin(1 - 1e-10, pmax(1e-10, inv_Logistic_arg))
    x <- 0.5 + sigma * qlogis(inv_Logistic_arg)
  } else {
    # Right reflected truncated exponential component
    u_rel <- (theta_q - w2 - w1) / w3
    exp_neg_lambda <- exp(-lambda)
    Z_E <- 1 - exp_neg_lambda
    exp_arg <- u_rel * Z_E + exp_neg_lambda
    exp_arg <- pmax(1e-10, exp_arg)
    x <- 1 + (1 / lambda) * log(exp_arg)
  }
  
  return(x)
}

# Compute the inverse CDF for a batch of parameters
compute_inverse_cdf_batch <- function(batch) {
  batch$theta <- pmap_dbl(batch, ~ compute_inverse_cdf(..4, ..1, ..2, ..3)) # ..4 = theta_q, ..1 = w1, etc.
  return(batch)
}

generate_inverse_cdf_grid <- library(parallel)

# Function to compute the inverse CDF
compute_inverse_cdf <- function(theta_q, w1, lambda, sigma) {
  w2 <- w3 <- (1 - w1) / 2
  
  if (theta_q < w2) {
    u <- theta_q / w2
    Z_E <- 1 - exp(-lambda)
    x <- - (1 / lambda) * log(1 - u * Z_E)
  } else if (theta_q < (w2 + w1)) {
    u <- (theta_q - w2) / w1
    Logistic_L <- plogis((0 - 0.5) / sigma)
    Logistic_U <- plogis((1 - 0.5) / sigma)
    inv_Logistic_arg <- u * (Logistic_U - Logistic_L) + Logistic_L
    inv_Logistic_arg <- pmin(1 - 1e-10, pmax(1e-10, inv_Logistic_arg))
    x <- 0.5 + sigma * qlogis(inv_Logistic_arg)
  } else {
    u_rel <- (theta_q - w2 - w1) / w3
    exp_neg_lambda <- exp(-lambda)
    Z_E <- 1 - exp_neg_lambda
    exp_arg <- u_rel * Z_E + exp_neg_lambda
    exp_arg <- pmax(1e-10, exp_arg)
    x <- 1 + (1 / lambda) * log(exp_arg)
  }
  
  return(x)
}


library(parallel)

# Function to compute the inverse CDF
compute_inverse_cdf <- function(theta_q, w1, lambda, sigma) {
  w2 <- w3 <- (1 - w1) / 2
  
  if (theta_q < w2) {
    u <- theta_q / w2
    Z_E <- 1 - exp(-lambda)
    x <- - (1 / lambda) * log(1 - u * Z_E)
  } else if (theta_q < (w2 + w1)) {
    u <- (theta_q - w2) / w1
    Logistic_L <- plogis((0 - 0.5) / sigma)
    Logistic_U <- plogis((1 - 0.5) / sigma)
    inv_Logistic_arg <- u * (Logistic_U - Logistic_L) + Logistic_L
    inv_Logistic_arg <- pmin(1 - 1e-10, pmax(1e-10, inv_Logistic_arg))
    x <- 0.5 + sigma * qlogis(inv_Logistic_arg)
  } else {
    u_rel <- (theta_q - w2 - w1) / w3
    exp_neg_lambda <- exp(-lambda)
    Z_E <- 1 - exp_neg_lambda
    exp_arg <- u_rel * Z_E + exp_neg_lambda
    exp_arg <- pmax(1e-10, exp_arg)
    x <- 1 + (1 / lambda) * log(exp_arg)
  }
  
  return(x)
}

# Generate the 5D array
generate_inverse_cdf_grid <- function(
    n_w1 = 21, w1_min = 0, w1_max = 1,
    n_lambda = 20, lambda_min_log = log(10), lambda_max_log = log(100),
    n_sigma = 20, sigma_min_log = log(0.01), sigma_max_log = log(0.5),
    n_theta_q = 50, theta_q_min = 0, theta_q_max = 1,
    n_cores = detectCores() - 1
) {
  # Create parameter grids
  w1_grid <- seq(w1_min, w1_max, length.out = n_w1)
  lambda_grid <- exp(seq(lambda_min_log, lambda_max_log, length.out = n_lambda))
  sigma_grid <- exp(seq(sigma_min_log, sigma_max_log, length.out = n_sigma))
  theta_q_grid <- seq(theta_q_min, theta_q_max, length.out = n_theta_q)
  
  # Initialize a 5D array
  inverse_cdf_array <- array(
    NA, dim = c(n_w1, n_lambda, n_sigma, n_theta_q),
    dimnames = list(
      w1 = w1_grid,
      lambda = lambda_grid,
      sigma = sigma_grid,
      theta_q = theta_q_grid
    )
  )
  
  # Generate all parameter combinations for w1, lambda, and sigma
  parameter_indices <- expand.grid(
    i_w1 = seq_len(n_w1),
    i_lambda = seq_len(n_lambda),
    i_sigma = seq_len(n_sigma)
  )
  
  # Function to compute a slice for a single combination of w1, lambda, sigma
  compute_slice <- function(indices) {
    i_w1 <- indices$i_w1
    i_lambda <- indices$i_lambda
    i_sigma <- indices$i_sigma
    w1 <- w1_grid[i_w1]
    lambda <- lambda_grid[i_lambda]
    sigma <- sigma_grid[i_sigma]
    
    # Compute theta for all theta_q values
    theta_values <- sapply(theta_q_grid, compute_inverse_cdf, w1 = w1, lambda = lambda, sigma = sigma)
    list(indices = indices, theta_values = theta_values)
  }
  
  # Compute in parallel
  results <- mclapply(seq_len(nrow(parameter_indices)), function(row) {
    compute_slice(parameter_indices[row, ])
  }, mc.cores = n_cores)
  
  # Fill the 5D array with results
  for (res in results) {
    i_w1 <- res$indices$i_w1
    i_lambda <- res$indices$i_lambda
    i_sigma <- res$indices$i_sigma
    inverse_cdf_array[i_w1, i_lambda, i_sigma, ] <- res$theta_values
  }
  
  return(inverse_cdf_array)
}

# Generate the 5D array
generate_inverse_cdf_grid_array <- function(
  n_w1 = 21, w1_min = 0, w1_max = 1,
  n_lambda = 20, lambda_min_log = log(10), lambda_max_log = log(100),
  n_sigma = 20, sigma_min_log = log(0.01), sigma_max_log = log(0.5),
  n_theta_q = 50, theta_q_min = 0, theta_q_max = 1,
  n_cores = detectCores() - 1
) {
  # Create parameter grids
  w1_grid <- seq(w1_min, w1_max, length.out = n_w1)
  lambda_grid <- exp(seq(lambda_min_log, lambda_max_log, length.out = n_lambda))
  sigma_grid <- exp(seq(sigma_min_log, sigma_max_log, length.out = n_sigma))
  theta_q_grid <- seq(theta_q_min, theta_q_max, length.out = n_theta_q)
  
  # Initialize a 5D array
  inverse_cdf_array <- array(
    NA, dim = c(n_w1, n_lambda, n_sigma, n_theta_q),
    dimnames = list(
      w1 = w1_grid,
      lambda = lambda_grid,
      sigma = sigma_grid,
      theta_q = theta_q_grid
    )
  )
  
  # Generate all parameter combinations for w1, lambda, and sigma
  parameter_indices <- expand.grid(
    i_w1 = seq_len(n_w1),
    i_lambda = seq_len(n_lambda),
    i_sigma = seq_len(n_sigma)
  )
  
  # Function to compute a slice for a single combination of w1, lambda, sigma
  compute_slice <- function(indices) {
    i_w1 <- indices$i_w1
    i_lambda <- indices$i_lambda
    i_sigma <- indices$i_sigma
    w1 <- w1_grid[i_w1]
    lambda <- lambda_grid[i_lambda]
    sigma <- sigma_grid[i_sigma]
    
    # Compute theta for all theta_q values
    theta_values <- sapply(theta_q_grid, compute_inverse_cdf, w1 = w1, lambda = lambda, sigma = sigma)
    list(indices = indices, theta_values = theta_values)
  }
  
  # Compute in parallel
  results <- mclapply(seq_len(nrow(parameter_indices)), function(row) {
    compute_slice(parameter_indices[row, ])
  }, mc.cores = n_cores)
  
  # Fill the 5D array with results
  for (res in results) {
    i_w1 <- res$indices$i_w1
    i_lambda <- res$indices$i_lambda
    i_sigma <- res$indices$i_sigma
    inverse_cdf_array[i_w1, i_lambda, i_sigma, ] <- res$theta_values
  }
  
  return(inverse_cdf_array)
}


# Generate the grid with default parameters
inverse_cdf_grid <- generate_inverse_cdf_grid(
  n_w1 = 50, w1_min = 0.1, w1_max = 0.99,
  n_lambda = 50, lambda_min_log = log(10), lambda_max_log = log(100),
  n_sigma = 50, sigma_min_log = log(0.01), sigma_max_log = log(0.5),
  n_theta_q = 100, theta_q_min = 0, theta_q_max = 1
)

# Save the grid to a file for later use
saveRDS(inverse_cdf_grid, file = "inverse_cdf_grid.rds")
