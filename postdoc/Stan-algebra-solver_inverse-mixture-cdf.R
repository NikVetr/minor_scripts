library(cmdstanr)
library(posterior)

model_index <- 2
model_paths <- c("~/scripts/minor_scripts/postdoc/Stan-algebra-solver_inverse-mixture-cdf.stan",
                "~/scripts/minor_scripts/postdoc/Stan-algebra-solver_inverse-mixture-logistic-center-cdf.stan")
model_path <- model_paths[model_index]
model <- cmdstan_model(model_path)

#### fit the model ####

# Define the data
mu <- 0.5       # Mean of the normal distribution
sigma <- 0.01    # scale parameter of the slab distribution
lambda <- 50   # Rate of the exponential distributions
w1 <- 0.5       # Weight of the normal component
wexp <- 1 - w1
w2 <- wexp / 2       # Weight of the exponential component
w3 <- wexp / 2       # Weight of the reflected exponential component
solver_index <- 2

# Ensure that the weights sum to 1
if (abs(w1 + w2 + w3 - 1) > 1e-6) {
  stop("Weights w1, w2, and w3 must sum to 1.")
}

# Create a list of data to pass to Stan
data <- list(
  mu = mu,
  sigma = sigma,
  lambda = lambda,
  w1 = w1,
  w2 = w2,
  w3 = w3,
  solver_index = solver_index
)

# Run MCMC sampling
niter <- 1E4
fit <- model$sample(
  data = data,
  seed = 123,               # For reproducibility
  chains = 4,               # Number of Markov chains
  parallel_chains = 4,      # Number of chains to run in parallel
  iter_warmup = niter / 4,        # Number of warm-up iterations per chain
  iter_sampling = niter * 3 / 4,     # Number of sampling iterations per chain
  refresh = 500             # Show progress every 100 iterations
)

# Extract the samples for theta
theta_samples_mcmc <- as_draws_df(fit$draws("theta"))$theta

# compare to results from direct simulation
n_samples <- niter

# Vectorized sampling of mixture components
component <- sample(1:3, size = n_samples, replace = TRUE, prob = c(w1, w2, w3))

# Initialize theta vector
theta_samples_R <- numeric(n_samples)

# Function to sample from a truncated normal distribution using inverse CDF
sample_truncnorm <- function(n, mu, sigma, lower, upper) {
  lower_p <- pnorm(lower, mean = mu, sd = sigma)  # CDF at lower bound
  upper_p <- pnorm(upper, mean = mu, sd = sigma)  # CDF at upper bound
  u <- runif(n, min = lower_p, max = upper_p)     # Rescaled uniform quantile
  qnorm(u, mean = mu, sd = sigma)                # Inverse CDF
}

sample_trunclogistic <- function(n, mu, sigma, lower, upper) {
  lower_p <- plogis(lower, location = mu, scale = sigma)  # CDF at lower bound
  upper_p <- plogis(upper, location = mu, scale = sigma)  # CDF at upper bound
  u <- runif(n, min = lower_p, max = upper_p)     # Rescaled uniform quantile
  qlogis(u, location = mu, scale = sigma)                # Inverse CDF
}

# Function to sample from a truncated exponential distribution using inverse CDF
sample_truncexp <- function(n, rate, lower, upper) {
  lower_p <- pexp(lower, rate = rate)  # CDF at lower bound
  upper_p <- pexp(upper, rate = rate)  # CDF at upper bound
  u <- runif(n, min = lower_p, max = upper_p)  # Rescaled uniform quantile
  qexp(u, rate = rate)                        # Inverse CDF
}

# Component 1: Truncated Normal
idx1 <- which(component == 1)
n1 <- length(idx1)
if (n1 > 0) {
  if(grepl("logistic", model_path)){
    theta_samples_R[idx1] <- sample_trunclogistic(n1, mu, sigma, lower = 0, upper = 1)
  } else {
    theta_samples_R[idx1] <- sample_truncnorm(n1, mu, sigma, lower = 0, upper = 1)
  }
  
}

# Component 2: Truncated Exponential
idx2 <- which(component == 2)
n2 <- length(idx2)
if (n2 > 0) {
  theta_samples_R[idx2] <- sample_truncexp(n2, rate = lambda, lower = 0, upper = 1)
}

# Component 3: Reflected Truncated Exponential
idx3 <- which(component == 3)
n3 <- length(idx3)
if (n3 > 0) {
  theta_samples_R[idx3] <- 1 - sample_truncexp(n3, rate = lambda, lower = 0, upper = 1)
}


# Plot the posterior distribution of theta
hist(theta_samples_R, breaks = 0:100/100, probability = TRUE, 
     main = "Simulated Posterior Distribution of theta",
     xlab = expression(theta), col = "skyblue", border = "black")
hist(theta_samples_mcmc, breaks = 0:100/100, probability = TRUE, add = T,
     xlab = expression(theta), col = adjustcolor("red", 0.5), border = "black")
legend(x = 0.8, y = par("usr")[3] - diff(par("usr")[3:4])/3, 
       col = c("skyblue", "red"), pch = 15, legend = c("direct", "mcmc"),
       xpd = NA, bty = "n")
