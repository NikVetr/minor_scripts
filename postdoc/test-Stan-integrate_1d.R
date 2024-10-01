# Load cmdstanr package
library(cmdstanr)

# Probability Density Function (PDF) for the Generalized Laplace Distribution with sigma
dlap <- function(x, mu, sigma) {
  lambda <- sigma / sqrt(2)
  return((1 / (2 * lambda)) * exp(-abs(x - mu) / lambda))
}

# Cumulative Distribution Function (CDF) for the Generalized Laplace Distribution with sigma
plap <- function(x, mu, sigma) {
  lambda <- sigma / sqrt(2)
  sapply(x, function(xi){
    if (xi < mu) {
      return(0.5 * exp((xi - mu) / lambda))
    } else {
      return(1 - 0.5 * exp(-(xi - mu) / lambda))
    }  
  })
}

# Random Generation for the Generalized Laplace Distribution with sigma
rlap <- function(n, mu, sigma) {
  lambda <- sigma / sqrt(2)
  u <- runif(n) - 0.5  # Uniform random numbers centered at 0
  return(mu - lambda * sign(u) * log(1 - 2 * abs(u)))
}

# Quantile Function for the Generalized Laplace Distribution with sigma
qlap <- function(p, mu, sigma) {
  lambda <- sigma / sqrt(2)
  sapply(p, function(pi){
    if (pi < 0.5) {
      return(mu + lambda * log(2 * pi))
    } else {
      return(mu - lambda * log(2 * (1 - pi)))
    }  
  })
}



#### just a single integral ####

#write basic model
stan_code_1 <- "
functions {
  real f(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    return exp(-square(x));  // define f(x) = exp(-x^2)
  }
}

data {
  real lower_bound;
  real upper_bound;
}

generated quantities {
  real integral;
  integral = integrate_1d(f, lower_bound, upper_bound, rep_array(0.0, 0), rep_array(0.0, 0), rep_array(0, 0));
}
"

# Compile the Stan model
mod_1 <- cmdstan_model(write_stan_file(stan_code_1))

# Define data for the model
dat_1 <- list(
  lower_bound = 0.0,
  upper_bound = 1.0
)

# Run the model
fit_1 <- mod_1$sample(data = dat_1, chains = 1, iter_sampling = 1, fixed_param = TRUE)

# Extract the result
samps_1 <- fit_1$draws(format = "df")

# Show the integral result
print(samps_1)

#validate externally

# Define the function to integrate in R
f_r <- function(x) {
  exp(-x^2)
}

# Perform numerical integration using R's integrate function
result_r <- integrate(f_r, 0, 1)

# Show the R integral result and compare
print(paste("Stan result:", samps_1$integral))
print(paste("R result:", result_r$value))


#### sample from a target distribution ####

#integration with asinh-tanh quadruture
stan_code_2 <- "
functions {
  real normal_integral(real sigma, real xc, array[] real theta, data array[] real x_r, data array[] int x_i) {
    real y = theta[1];
    real mu = x_r[1];
    real tau = x_r[2];  // Scale parameter for the half-normal prior
    
    // Ensure sigma is positive
    if (sigma <= 0) {
      return 0.0;  // The half-normal is defined for sigma > 0
    }
    
    // Normal likelihood for y given mu and sigma
    real normal_log_density = normal_lpdf(y | mu, sigma);
    
    // Half-normal prior on sigma (log density)
    real half_normal_log_density = normal_lpdf(sigma | 0, tau) + log(2);  // Half-normal adjustment
    
    // Return the exponentiated sum of log densities
    return exp(normal_log_density + half_normal_log_density);
  }
}

data {
  real mu;     // Mean of the normal distribution
  real tau;    // Scale parameter of the half-normal prior on sigma
  real lb;    // lower bound of integration
  real ub;    // upper bound of integration
}

parameters {
  real y;
}

model {
  // Perform the integration
  target += log(integrate_1d(normal_integral, lb, ub, {y}, {mu, tau}, {0}, 1E-8));
}
"

#integration with parameter expansion / data augmentation
stan_code_2b <- "
data {
  real mu;     // Mean of the normal distribution
  real tau;    // Scale parameter of the half-normal prior on sigma
  real lb;    // lower bound of integration
  real ub;    // upper bound of integration
}

parameters {
  real<lower=lb, upper=ub> sigma;
  real y;
}

model {
  sigma ~ normal(0, tau);
  y ~ normal(mu, sigma);
}
"

#integration with discrete approximation (riemann sum)
stan_code_2c <- "
data {
  real mu;              // Mean of the normal distribution
  real tau;             // Scale parameter of the half-normal prior on sigma
  real lb;              // Lower bound of integration
  real ub;              // Upper bound of integration
  int<lower=1> nblocks; // Number of blocks for numerical integration
}

transformed data {
  real delta_sigma = (ub - lb) / nblocks; // Width of each sigma interval
  vector[nblocks] sigma;                   // Array to store sigma values
  
  for (i in 1:nblocks) {
    // Midpoint of each sigma interval
    sigma[i] = lb + (i - 0.5) * delta_sigma;
  }
}

parameters {
  real y; // Observation to be sampled
}

model {
  vector[nblocks] log_p_integrand;
  
  for (i in 1:nblocks) {
    // Compute the log-likelihood of y given sigma[i]
    real log_p_y_given_sigma = normal_lpdf(y | mu, sigma[i]);
    // Compute the log of the half-normal prior at sigma[i]
    real log_p_sigma = normal_lpdf(sigma[i] | 0, tau) + log(2); // Half-normal prior
    // Sum the log probabilities
    log_p_integrand[i] = log_p_y_given_sigma + log_p_sigma;
  }
  
  // Compute the log of the marginal likelihood of y
  target += log(delta_sigma) + log_sum_exp(log_p_integrand);
}
"


#integration with discrete approximation (riemann sum) in transformed data
stan_code_2c.2 <- "
data {
  real mu;              // Mean of the normal distribution
  real tau;             // Scale parameter of the half-normal prior on sigma
  real lb;              // Lower bound of integration
  real ub;              // Upper bound of integration
  int<lower=1> nblocks; // Number of blocks for numerical integration
}

transformed data {
  //  real delta_sigma = (ub - lb) / nblocks; // Width of each sigma interval
  //  vector[nblocks] sigma;                   // Array to store sigma values
    
  //  for (i in 1:nblocks) {
  //    // Midpoint of each sigma interval
  //    sigma[i] = lb + (i - 0.5) * delta_sigma;
  //  }
  
  real cdf_lb = normal_cdf(lb | 0, tau); // CDF at lower bound
  real cdf_ub = normal_cdf(ub | 0, tau); // CDF at upper bound
  real delta_cdf = (cdf_ub - cdf_lb) / nblocks;  // Probability mass per block
  vector[nblocks] sigma;                         // Array to store sigma values

  for (i in 1:nblocks) {
    // Calculate the midpoint of each probability block in CDF space
    real cdf_mid = cdf_lb + (i - 0.5) * delta_cdf;
    
    // Compute sigma from the inverse CDF (inverse of half-normal)
    sigma[i] = tau * inv_Phi(cdf_mid);
  }
}

parameters {
  real y; // Observation to be sampled
}

model {
  vector[nblocks] log_p_integrand;
  
  for (i in 1:nblocks) {
    // Compute the log-likelihood of y given sigma[i]
    real log_p_y_given_sigma = normal_lpdf(y | mu, sigma[i]);
    // Compute the log of the half-normal prior at sigma[i]
    // real log_p_sigma = normal_lpdf(sigma[i] | 0, tau) + log(2); // Half-normal prior
    // Sum the log probabilities
    log_p_integrand[i] = log_p_y_given_sigma; // + log_p_sigma;
  }
  
  // Compute the log of the marginal likelihood of y
  target += log_sum_exp(log_p_integrand);
}
"


#integration with discrete approximation (riemann sum) in function
stan_code_2d <- "
functions {
  real discrete_approx_normal_mixture_logpdf(
      real y, real mu, real tau, real lb, real ub, int nblocks) {
    // Width of each sigma interval
    real delta_sigma = (ub - lb) / nblocks;
    vector[nblocks] log_p_integrand;
    real cdf_lb = normal_lcdf(lb | 0, tau); // CDF at lower bound
    real cdf_ub = normal_lcdf(ub | 0, tau); // CDF at upper bound
    real trunc_log_p = log_diff_exp(cdf_ub, cdf_lb);

    // Compute sigma and the corresponding log densities
    for (i in 1:nblocks) {
      real sigma = lb + (i - 0.5) * delta_sigma; // Midpoint of each sigma interval
      // Compute the log-likelihood of y given sigma
      real log_p_y_given_sigma = normal_lpdf(y | mu, sigma);
      // Compute the log of the half-normal prior at sigma
      real log_p_sigma = normal_lpdf(sigma | 0, tau) + log(2); // Half-normal prior
      // Sum the log probabilities
      log_p_integrand[i] = log_p_y_given_sigma + log_p_sigma + trunc_log_p;
    }
    
    // Compute the log of the marginal likelihood of y
    return log(delta_sigma) + log_sum_exp(log_p_integrand) ;
  }
}

data {
  real mu;              // Mean of the normal distribution
  real tau;             // Scale parameter of the half-normal prior on sigma
  real lb;              // Lower bound of integration
  real ub;              // Upper bound of integration
  int<lower=1> nblocks; // Number of blocks for numerical integration
}

parameters {
  real y; // Observation to be sampled
}

model {
  // Call the discrete approximation function
  target += discrete_approx_normal_mixture_logpdf(y, mu, tau, lb, ub, nblocks);
}
"


# Compile the Stan model
mod_2 <- cmdstan_model(write_stan_file(stan_code_2))
mod_2b <- cmdstan_model(write_stan_file(stan_code_2b))
mod_2c <- cmdstan_model(write_stan_file(stan_code_2c))
mod_2d <- cmdstan_model(write_stan_file(stan_code_2d))

# Define data for the model
dat_2 <- list(
  mu = 1.3,        # location
  tau = 0.6,     # scale of scale
  lb = 0.7,
  ub = 1.9,
  nblocks = 200
)

niter <- 1E4

# Run the model
fit_2 <- mod_2$sample(data = dat_2, chains = 4, parallel_chains = 4, refresh = 1E3,
                      iter_sampling = niter, iter_warmup = niter, thin = 10)
summ_2 <- fit_2$summary()
summ_2

fit_2b <- mod_2b$sample(data = dat_2, chains = 4, parallel_chains = 4, refresh = 1E3,
                        iter_sampling = niter, iter_warmup = niter, thin = 10)
summ_2b <- fit_2b$summary()
summ_2b

fit_2c <- mod_2c$sample(data = dat_2, chains = 4, parallel_chains = 4, refresh = 1E3,
                        iter_sampling = niter, iter_warmup = niter, thin = 10)
summ_2c <- fit_2c$summary()
summ_2c

fit_2d <- mod_2d$sample(data = dat_2, chains = 4, parallel_chains = 4, refresh = 1E3,
                        iter_sampling = niter, iter_warmup = niter, thin = 10)
summ_2d <- fit_2d$summary()
summ_2d

# Extract the result
samps_2 <- fit_2$draws(format = "df")
samps_2b <- fit_2b$draws(format = "df")
samps_2c <- fit_2c$draws(format = "df")
samps_2d <- fit_2d$draws(format = "df")

#compare to direct sampling in R
n_samples <- niter  # Number of samples to generate

# Simulate sigma from half-normal distribution
uniform_rv <- runif(n_samples, 0, 1)
bounds_quantiles <- pnorm(c(dat_2$lb, dat_2$ub), sd = dat_2$tau)
normal_quantiles <- bounds_quantiles[1] + diff(bounds_quantiles) * uniform_rv
sigma_samples <- qnorm(normal_quantiles, sd = dat_2$tau)

# Simulate y from normal distribution using the sigma samples
y_samples_r <- rnorm(n_samples, mean = dat_2$mu, sd = sigma_samples)

#compare Stan samples to direct samples
mean(y_samples_r)
mean(samps_2$y)
mean(samps_2b$y)
mean(samps_2c$y)
mean(samps_2d$y)

var(y_samples_r)
var(samps_2$y)
var(samps_2b$y)
var(samps_2c$y)
var(samps_2d$y)

breaks = -20:20/2
hist(y_samples_r, add = F, col = adjustcolor(3,0.2), freq = F, breaks = breaks)
hist(samps_2$y, freq = F, breaks = breaks, col = adjustcolor(3,0.2), add = T)
hist(samps_2b$y, add = T, col = adjustcolor(4,0.2), freq = F, breaks = breaks)
hist(samps_2c$y, add = T, col = adjustcolor(5,0.2), freq = F, breaks = breaks)
hist(samps_2d$y, add = T, col = adjustcolor(6,0.2), freq = F, breaks = breaks)

#### sample from trunc laplace binomial ####

#use discrete approximation
stan_code_3 <- "
functions {

  real discrete_approx_laplace_binomial_logpmf(
      array[] int count, array[] int total, real mu, real sigma, real lb, real ub, int nblocks) {
    
    // Width of each sigma interval
    int n = size(count);
    real delta_theta = (ub - lb) / nblocks;
    vector[nblocks] theta;
    vector[nblocks] log_p_theta;
    vector[n] integrated_lp;
    
    // truncated distribution normalization factor
    real cdf_lb = double_exponential_lcdf(lb | mu, sigma); // CDF at lower bound
    real cdf_ub = double_exponential_lcdf(ub | mu, sigma); // CDF at upper bound
    real trunc_log_p = log_diff_exp(cdf_ub, cdf_lb);
    
    // Compute theta and the corresponding log densities
    for (i in 1:nblocks) {
      theta[i] = lb + (i - 0.5) * delta_theta; // Midpoint of each theta interval
      log_p_theta[i] = double_exponential_lpdf(theta[i] | mu, sigma) - trunc_log_p;
    }
    real lse_theta = log_sum_exp(log_p_theta);
    
    // Compute the log-likelihood of data given theta
    for(j in 1:n){
      vector[nblocks] log_p_count_given_theta;
      for (i in 1:nblocks) {
        log_p_count_given_theta[i] = binomial_lpmf(count[j] | total[j], theta[i]);
      }
      integrated_lp[j] = log_sum_exp(log_p_theta + log_p_count_given_theta) - lse_theta;
    }
    
    // Compute the log of the marginal likelihood of count
    return sum(integrated_lp);
  }
}

data {
  real lb;                      // Lower bound of integration
  real ub;                      // Upper bound of integration
  int<lower=1> nblocks;         // Number of blocks for numerical integration
  int<lower=1> n;               // Number of observations
  array[n] int<lower=1> total;  // Number of trials for the binomial distribution
  array[n] int<lower=0> count;  // Number of successes
}

parameters {
  real mu; // Mean of the laplace distribution
  real<lower=0> sigma; // Scale parameter of laplace distribution
}

model {
  // priors
  mu ~ std_normal();
  sigma ~ std_normal();

  // likelihood
  target += discrete_approx_laplace_binomial_logpmf(count, total, mu, sigma, 
                                                    lb, ub, nblocks);
}
"

#work on raw scale?
stan_code_3a <- "
functions {

  real discrete_approx_laplace_binomial_logpmf(
      array[] int count, array[] int total, real mu, real sigma, real lb, real ub, int nblocks) {
    
    // Width of each sigma interval
    int n = size(count);
    real delta_theta = (ub - lb) / nblocks;
    vector[nblocks] theta;
    vector[nblocks] p_theta;
    vector[n] integrated_p;
    
    // truncated distribution normalization factor
    real cdf_lb = double_exponential_cdf(lb | mu, sigma); // CDF at lower bound
    real cdf_ub = double_exponential_cdf(ub | mu, sigma); // CDF at upper bound
    real trunc_p = cdf_ub - cdf_lb;
    
    // Compute theta and the corresponding log densities
    for (i in 1:nblocks) {
      theta[i] = lb + (i - 0.5) * delta_theta; // Midpoint of each theta interval
      p_theta[i] = exp(double_exponential_lpdf(theta[i] | mu, sigma)) / trunc_p;
    }
    real sum_theta_mass = sum(p_theta);
    
    // Compute the log-likelihood of data given theta
    for(j in 1:n){
      vector[nblocks] p_count_given_theta;
      for (i in 1:nblocks) {
        p_count_given_theta[i] = exp(binomial_lpmf(count[j] | total[j], theta[i]));
      }
      integrated_p[j] = sum(p_theta .* p_count_given_theta / sum_theta_mass);
    }
    
    // Compute the log of the marginal likelihood of count
    return log(prod(integrated_p));
  }
}

data {
  real lb;                      // Lower bound of integration
  real ub;                      // Upper bound of integration
  int<lower=1> nblocks;         // Number of blocks for numerical integration
  int<lower=1> n;               // Number of observations
  array[n] int<lower=1> total;  // Number of trials for the binomial distribution
  array[n] int<lower=0> count;  // Number of successes
}

parameters {
  real mu; // Mean of the laplace distribution
  real<lower=0> sigma; // Scale parameter of laplace distribution
}

model {
  // priors
  mu ~ std_normal();
  sigma ~ std_normal();

  // likelihood
  target += discrete_approx_laplace_binomial_logpmf(count, total, mu, sigma, 
                                                    lb, ub, nblocks);
}
"

#approximate through parameter expansion
stan_code_3b <- "
data {
  real lb;              // Lower bound of integration
  real ub;              // Upper bound of integration
  int<lower=1> nblocks; // Number of blocks for numerical integration
  int<lower=1> n;               // Number of observations
  array[n] int<lower=1> total;  // Number of trials for the binomial distribution
  array[n] int<lower=0> count;  // Number of successes
}

parameters {
  array[n] real<lower=lb, upper=ub> theta; // binomial probability
  real mu; // Mean of the laplace distribution
  real<lower=0> sigma; // Scale parameter of laplace distribution
}

model {
  // priors
  mu ~ std_normal();
  sigma ~ std_normal();
  theta ~ double_exponential(mu, sigma);

  // likelihood
  count ~ binomial(total, theta);
}
"

stan_code_3c <- "
functions {
  
  real laplace_binomial_integral(real x, real xc, array[] real theta, data array[] real x_r, data array[] int x_i) {
    real mu = theta[1];
    real sigma = theta[2];
    int count = x_i[1];
    int total = x_i[2];
    
    real binomial_log_mass = binomial_lpmf(count | total, x);
    real laplace_log_density = double_exponential_lpdf(x | mu, sigma);
    
    // Return the exponentiated sum of log densities
    return exp(binomial_log_mass + laplace_log_density);
  }
}

data {
  real lb;              // Lower bound of integration
  real ub;              // Upper bound of integration
  int<lower=1> nblocks; // Number of blocks for numerical integration
  int<lower=1> n;               // Number of observations
  array[n] int<lower=1> total;  // Number of trials for the binomial distribution
  array[n] int<lower=0> count;  // Number of successes
}

parameters {
  real mu; // Mean of the laplace distribution
  real<lower=0> sigma; // Scale parameter of laplace distribution
}

model {
  // Perform the integration
  for(i in 1:n){
    target += log(integrate_1d(laplace_binomial_integral, lb, ub, {mu, sigma}, {1.0}, 
                  {count[i], total[i]}, 1E-8));
  }
}
"



# Compile the Stan model
mod_3 <- cmdstan_model(write_stan_file(stan_code_3))
mod_3a <- cmdstan_model(write_stan_file(stan_code_3a))
mod_3b <- cmdstan_model(write_stan_file(stan_code_3b))
mod_3c <- cmdstan_model(write_stan_file(stan_code_3c))

# Generate some sample data
n <- 20
mu <- 0.5
sigma <- 0.1
unif_rv <- runif(n, 0, 1)
bounds <- c(lb = 0.4, ub = 0.8)
bounds_quantiles <- plap(x = bounds, sigma = sigma, mu = mu)
sample_quantiles <- bounds_quantiles[1] + unif_rv * diff(bounds_quantiles)
theta <- qlap(p = sample_quantiles, mu = mu, sigma = sigma)
total <- 100
count <- rbinom(n = n, size = total, prob = theta)

#put into a list for Stan
dat_3 <- list(
  n = n,
  count = count,
  total = rep(total, n),
  nblocks = 200,      # Number of discrete points
  lb = bounds["lb"],           # Lower bound for theta
  ub = bounds["ub"]            # Upper bound for theta
)

niter <- 1E3

# Run the model
fit_3 <- mod_3$sample(data = dat_3, chains = 4, parallel_chains = 4, refresh = ceiling(niter / 20),
                      iter_sampling = niter, iter_warmup = niter, thin = ceiling(niter / 1E3))
summ_3 <- fit_3$summary()
summ_3

fit_3a <- mod_3a$sample(data = dat_3, chains = 4, parallel_chains = 4, refresh = ceiling(niter / 20),
                        iter_sampling = niter, iter_warmup = niter, thin = ceiling(niter / 1E3))
summ_3a <- fit_3a$summary()
summ_3a

fit_3c <- mod_3c$sample(data = dat_3, chains = 4, parallel_chains = 4, refresh = ceiling(niter / 20),
                        iter_sampling = niter, iter_warmup = niter, thin = ceiling(niter / 1E3))
summ_3c <- fit_3c$summary()
summ_3c

fit_3b <- mod_3b$sample(data = dat_3, chains = 4, parallel_chains = 4, refresh = ceiling(niter / 20),
                      iter_sampling = niter, iter_warmup = niter, thin = ceiling(niter / 1E3))
# summ_3b <- fit_3b$summary()
# summ_3b


#example model output
samps_3 <- fit_3$draws(format = "df")
samps_3a <- fit_3a$draws(format = "df")
samps_3b <- fit_3b$draws(format = "df")
nbins <- 50

range_mu <- range(samps_3$mu, samps_3b$mu)
breaks_mu <- seq(range_mu[1], range_mu[2], length.out = nbins+1)
hist(samps_3$mu, freq = F, col = adjustcolor(2,0.2), breaks = breaks_mu)
hist(samps_3a$mu, add = T, col = adjustcolor(4,0.2), freq = F, breaks = breaks_mu)
hist(samps_3b$mu, add = T, col = adjustcolor(3,0.2), freq = F, breaks = breaks_mu)
abline(v = mu, col = 4, lwd = 2)

range_sigma <- range(samps_3$sigma, samps_3b$sigma)
breaks_sigma <- seq(range_sigma[1], range_sigma[2], length.out = nbins+1)
hist(samps_3$sigma, freq = F, col = adjustcolor(2,0.2), breaks = breaks_sigma)
hist(samps_3a$sigma, freq = F, col = adjustcolor(4,0.2), breaks = breaks_sigma, add = T)
hist(samps_3b$sigma, add = T, col = adjustcolor(3,0.2), freq = F, breaks = breaks_sigma)
abline(v = sigma, col = 4, lwd = 2)

