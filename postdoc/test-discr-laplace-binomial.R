# Load the cmdstanr library
library(cmdstanr)

#functions for simulation
# Probability Density Function (PDF) for the Generalized Laplace Distribution
dlap <- function(x, mu, lambda) {
  return((1 / (2 * lambda)) * exp(-abs(x - mu) / lambda))
}

# Cumulative Distribution Function (CDF) for the Generalized Laplace Distribution
plap <- function(x, mu, lambda) {
  sapply(x, function(xi){
    if (xi < mu) {
      return(0.5 * exp((xi - mu) / lambda))
    } else {
      return(1 - 0.5 * exp(-(xi - mu) / lambda))
    }  
  })
}

# Random Generation for the Generalized Laplace Distribution
rlap <- function(n, mu, lambda) {
  u <- runif(n) - 0.5  # Uniform random numbers centered at 0
  return(mu - lambda * sign(u) * log(1 - 2 * abs(u)))
}

# Quantile Function for the Generalized Laplace Distribution
qlap <- function(p, mu, lambda) {
  sapply(p, function(pi){
    if (pi < 0.5) {
      return(mu + lambda * log(2 * pi))
    } else {
      return(mu - lambda * log(2 * (1 - pi)))
    }  
  })
}


# Define the Stan model
stan_code_discr <- "
functions {
  // Log probability mass function for a binomial averaged over a Laplace distribution
  real discr_laplace_binomial_logpmf(int count, int total, real mu, real lambda, int N_pts, real lb, real ub) {
    
    //initialize variables
    real delta = (ub - lb) / (N_pts - 1); // Step size for equally spaced points
    vector[N_pts] theta;
    vector[N_pts] log_weights;
    vector[N_pts] log_probs;
    
    // Compute log probabilities and adjusted log weights
    for (i in 1:N_pts) {
      theta[i] = lb + (i - 1) * delta;
      log_weights[i] = double_exponential_lpdf(theta[i] | mu, lambda) - log_normalization;
      log_probs[i] = binomial_lpmf(count | total, theta[i]);
    }
    
    // Compute log-sum-exp for weighted log-probabilities
    real log_weighted_sum = log_sum_exp(log_probs + log_weights);  // log(p_i) + log(w_i) for each valid point
    real log_weights_sum = log_sum_exp(log_weights);  // Sum of log-weights

    return log_weighted_sum - log_weights_sum;  // Return the normalized log probability
  }
}

data {
  int<lower=1> n;
  array[n] int<lower=0> total;           // Number of trials for the binomial distribution
  array[n] int<lower=0> count;           // Number of successes
  int<lower=1> N_pts;       // Number of points for the discrete Laplace approximation
  real lb;                  // Lower bound for the probability parameter
  real ub;                  // Upper bound for the probability parameter
}

parameters {
  real mu;
  real<lower=0> lambda;
}

model {
  // Priors for mu and lambda
  mu ~ normal(0.5, 0.5);
  lambda ~ exponential(10);

  // Discrete approximation over Laplace-distributed probabilities
  for(i in 1:n){
    target += discr_laplace_binomial_logpmf(count[i], 
                total[i], mu, lambda, N_pts, lb, ub);
  }
}
"

stan_code_direct <- "
data {
  int<lower=1> n;
  array[n] int<lower=0> total;           // Number of trials for the binomial distribution
  array[n] int<lower=0> count;           // Number of successes
  real lb;                  // Lower bound for the probability parameter
  real ub;                  // Upper bound for the probability parameter
}

parameters {
  vector<lower=lb, upper=ub>[n] theta;
  real mu;
  real<lower=0> lambda;
}

model {
  // Priors for mu and lambda
  mu ~ normal(0.5, 0.5);
  lambda ~ exponential(10);
  
  //prior for theta
  theta ~ double_exponential(mu, lambda);

  // Binomial likelihood with direct Laplace-sampled theta
  count ~ binomial(total, theta);
}
"


# Compile the Stan model
mod_discr <- cmdstan_model(write_stan_file(stan_code_discr))
mod_direct <- cmdstan_model(write_stan_file(stan_code_direct))

# Generate some sample data
n <- 40
mu <- 0.6
lambda <- 0.1
unif_rv <- runif(n, 0, 1)
bounds <- c(lb = 0.5, ub = 0.8)
bounds_quantiles <- plap(x = bounds, lambda = lambda, mu = mu)
sample_quantiles <- bounds_quantiles[1] + unif_rv * diff(bounds_quantiles)
theta <- qlap(p = sample_quantiles, mu = mu, lambda = lambda)
total <- 100
count <- rbinom(n = n, size = total, prob = theta)

#put into a list for Stan
dat <- list(
  n = n,
  count = count,
  total = rep(total, n),
  N_pts = 5,      # Number of discrete points
  lb = bounds["lb"],           # Lower bound for theta
  ub = bounds["ub"]            # Upper bound for theta
)

# Fit the model
fit_discr <- mod_discr$sample(data = dat, iter_sampling = 1000, 
                              iter_warmup = 500, chains = 4, parallel_chains = 4)
fit_direct <- mod_direct$sample(data = dat, iter_sampling = 1000, 
                                iter_warmup = 500, chains = 4, parallel_chains = 4)

# Print the results
fit_discr$print()
fit_direct$print()

# Extract the samples
samps_discr <- data.frame(as_draws_df(fit_discr$draws()))
samps_direct <- data.frame(as_draws_df(fit_direct$draws()))

#plot results
mu_range <- range(c(samps_discr$mu, samps_direct$mu))
mu_breaks <- seq(mu_range[1], mu_range[2], length.out = 50)
hist_discr_mu <- hist(samps_discr$mu, plot = F, breaks = mu_breaks)
hist_direct_mu <- hist(samps_direct$mu, plot = F, breaks = mu_breaks)
plot(hist_discr_mu, col = adjustcolor(2, 0.5), 
     ylim = range(c(hist_discr_mu$counts, hist_direct_mu$counts)))
plot(hist_direct_mu, col = adjustcolor(3,0.5), add = T)
abline(v = mu, col = 4, lwd = 2, lty = 2)

lambda_range <- range(c(samps_discr$lambda, samps_direct$lambda))
lambda_breaks <- seq(lambda_range[1], lambda_range[2], length.out = 50)
hist_discr_lambda <- hist(samps_discr$lambda, plot = F, breaks = lambda_breaks)
hist_direct_lambda <- hist(samps_direct$lambda, plot = F, breaks = lambda_breaks)
plot(hist_discr_lambda, col = adjustcolor(2, 0.5), 
     ylim = range(c(hist_discr_lambda$counts, hist_direct_lambda$counts)))
plot(hist_direct_lambda, col = adjustcolor(3,0.5), add = T)
abline(v = lambda, col = 4, lwd = 2, lty = 2)
