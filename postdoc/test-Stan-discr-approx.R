# Load the cmdstanr library
library(cmdstanr)

# Define the Stan model
stan_code_discr <- "
functions {
  real discr_exp_normal(real x, real mu, real sigma, int N_pts, real lb, real ub) {
    
    //initialize variables
    real delta = (ub - lb) / (N_pts - 1); // Step size for equally spaced points
    vector[N_pts] theta;
    vector[N_pts] log_weights;
    vector[N_pts] log_probs;
    
    // Compute log probabilities and adjusted log weights
    for (i in 1:N_pts) {
      theta[i] = lb + (i - 1) * delta;
      //log_weights[i] = exponential_lpdf(theta[i] | 1 / (2 * sigma^2));
      log_weights[i] = exponential_lpdf(theta[i] | 1 / sigma);
      log_probs[i] = normal_lpdf(x | mu, theta[i]);
    }
    
    // Compute log-sum-exp for weighted log-probabilities
    real log_weighted_sum = log_sum_exp(log_probs + log_weights);  // log(p_i) + log(w_i) for each valid point
    real log_weights_sum = log_sum_exp(log_weights);  // Sum of log-weights

    return log_weighted_sum - log_weights_sum;  // Return the normalized log probability
  }
}

data{
  real mu;
  real<lower=0> sigma;
  int<lower=1> N_pts;
  real lb;
  real ub;
}

parameters {
  real x;
}

model {
  // Discrete approximation over exponentially distributed scale mixture of normals
  target += discr_exp_normal(x, mu, sigma, N_pts, lb, ub);
}
"

stan_code_direct <- "
data{
  real mu;
  real sigma;
  real lb;
  real ub;
}
parameters {
  real x;
}
model {
  // direct sampling of laplace
  x ~ double_exponential(mu, sigma);
}
"


# Compile the Stan model
mod_discr <- cmdstan_model(write_stan_file(stan_code_discr))
mod_direct <- cmdstan_model(write_stan_file(stan_code_direct))

#put into a list for Stan
dat <- list(
  mu = 3,
  sigma = 2,
  N_pts = 1000,
  lb = 0.001,
  ub = 10
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
x_range <- range(c(samps_discr$x, samps_direct$x))
x_breaks <- seq(x_range[1], x_range[2], length.out = 50)
hist_discr_x <- hist(samps_discr$x, plot = F, breaks = x_breaks)
hist_direct_x <- hist(samps_direct$x, plot = F, breaks = x_breaks)
plot(hist_discr_x, col = adjustcolor(2, 0.5), 
     ylim = range(c(hist_discr_x$counts, hist_direct_x$counts)))
plot(hist_direct_x, col = adjustcolor(3,0.5), add = T)




#try quadruture method?
library(statmod)

# Generate Gaussian-Legendre quadrature points and weights
N_pts <- 50  # Number of quadrature points (increase for more accuracy)
quad <- gauss.quad(N_pts, kind = "legendre")

# Transform the quadrature points to the correct interval for an exponential distribution
# Here we perform the change of variables necessary to adjust for the bounds of the exponential
quad_points <- -log(0.5 * (quad$nodes + 1))  # Map quadrature nodes to (0, infinity)
quad_weights <- quad$weights / quad_points   # Adjust weights for exponential distribution

# Data for Stan
dat_quad <- list(
  mu = 3,
  sigma = 2,
  N_pts = N_pts,
  quad_points = quad_points,
  quad_weights = quad_weights
)

# Define the Stan model
stan_code_quad <- "
functions {
  real discr_exp_normal_quadrature(real x, real mu, real sigma, int N_pts, vector quad_points, vector quad_weights) {
    
    vector[N_pts] log_weights;
    vector[N_pts] log_probs;
    
    // Compute log probabilities and adjusted log weights using quadrature points and weights
    for (i in 1:N_pts) {
      real theta = quad_points[i];  // Quadrature point for the standard deviation
      log_weights[i] = exponential_lpdf(theta | 1 / sigma); // Exponential prior on theta
      log_probs[i] = normal_lpdf(x | mu, theta);
    }
    
    // Compute log-sum-exp for weighted log-probabilities
    real log_weighted_sum = log_sum_exp(log_probs + log_weights + log(quad_weights));  // log(p_i) + log(w_i) + log(quad_weights[i])
    real log_weights_sum = log_sum_exp(log_weights + log(quad_weights));  // Sum of log-weights with quadrature weights

    return log_weighted_sum - log_weights_sum;  // Return the normalized log probability
  }
}

data {
  real mu;
  real<lower=0> sigma;
  int<lower=1> N_pts;
  vector[N_pts] quad_points;   // Quadrature points
  vector[N_pts] quad_weights;  // Quadrature weights
}

parameters {
  real x;
}

model {
  // Use the quadrature-based approximation for the mixture of normals
  target += discr_exp_normal_quadrature(x, mu, sigma, N_pts, quad_points, quad_weights);
}
"

# Compile and run the Stan model
mod_quad <- cmdstan_model(write_stan_file(stan_code_quad))
fit_quad <- mod_quad$sample(data = dat_quad, iter_sampling = 1000, iter_warmup = 500, chains = 4, parallel_chains = 4)

# Print the fit results
fit_quad$print()
fit_direct$print()
