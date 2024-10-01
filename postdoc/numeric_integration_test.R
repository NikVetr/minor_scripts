
#specify high level parameters
shapes <- c(a = 12, b = 17)
nsamps <- 10
count <- sample(1:200, nsamps, replace = T)
total <- sample(200:400, nsamps, replace = T)

#compute analytically
prod(extraDistr::dbbinom(x = count, size = total, 
                    alpha = shapes["a"], beta = shapes["b"]))

#compute numerically
n <- 100
theta <- seq(0,1,length.out=n)
beta_dens <- dbeta(x = theta, shape1 = shapes["a"], shape2 = shapes["b"])
beta_dens <- beta_dens / sum(beta_dens)
binom_prob <- do.call(rbind, lapply(1:n, function(i) dbinom(x = count, size = total, prob = theta[i])))
prod(apply(beta_dens * binom_prob, 2, sum))

#now do log scale

#compute numerically
beta_dens <- dbeta(x = theta, shape1 = shapes["a"], shape2 = shapes["b"], log = T) 
binom_prob <- do.call(rbind, lapply(1:n, function(i) dbinom(x = count, size = total, prob = theta[i], log = T)))

#this is effectively what the current function in stan is doing -- and it is wrong
# sum(matrixStats::logSumExp(beta_dens + apply(binom_prob, 1, sum)) - matrixStats::logSumExp(beta_dens))

#I need to do this instead
sum(
  sapply(1:nsamps, function(i){matrixStats::logSumExp(beta_dens + binom_prob[,i])}) - 
    matrixStats::logSumExp(beta_dens)
  )

#or
sum(matrixStats::colLogSumExps(binom_prob + beta_dens) - matrixStats::logSumExp(beta_dens))

#compute analytically
sum(extraDistr::dbbinom(x = count, size = total, 
                        alpha = shapes["a"], beta = shapes["b"], log = T))

#now get Stan to do it?
library(cmdstanr)

#write basic model
stan_code_1 <- "
functions {
  
  real discrete_approx_laplace_binomial_logpmf(
      array[] int count, array[] int total, real alpha, real beta, real lb, real ub, int nblocks) {
    
    // Width of each sigma interval
    int n = size(count);
    real delta_theta = (ub - lb) / nblocks;
    vector[nblocks] theta;
    vector[nblocks] log_p_theta;
    vector[n] integrated_lp;
    
    // truncated distribution normalization factor
    real cdf_lb = beta_lcdf(lb | alpha, beta); // CDF at lower bound
    real cdf_ub = beta_lcdf(ub | alpha, beta); // CDF at upper bound
    real trunc_log_p = log_diff_exp(cdf_ub, cdf_lb);
    
    // Compute theta and the corresponding log densities
    for (i in 1:nblocks) {
      theta[i] = lb + (i - 0.5) * delta_theta; // Midpoint of each theta interval
      log_p_theta[i] = beta_lpdf(theta[i] | alpha, beta) - trunc_log_p;
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
  real<lower=0> alpha;                      // Lower bound of integration
  real<lower=0> beta;                      // Upper bound of integration
  
}

generated quantities {
  real integral;
  integral = discrete_approx_laplace_binomial_logpmf(count, total, alpha, beta, 
                                                    lb, ub, nblocks);
}
"

# Compile the Stan model
mod_1 <- cmdstan_model(write_stan_file(stan_code_1))

# Define data for the model
bounds <- data.frame(lb = 0.3, ub = 0.6)
dat_1 <- list(
  n = nsamps,
  count = count,
  total = total,
  nblocks = 1000,
  lb = bounds$lb,
  ub = bounds$ub,
  alpha = shapes["a"],
  beta = shapes["b"]
)

# Run the model
fit_1 <- mod_1$sample(data = dat_1, chains = 1, iter_sampling = 1, fixed_param = TRUE)

# Extract the result
samps_1 <- fit_1$draws(format = "df")
samps_1$integral

#compute numerically
n <- 1000
theta <- seq(bounds$lb, bounds$ub, length.out=n)
beta_dens <- dbeta(x = theta, shape1 = shapes["a"], shape2 = shapes["b"], log = T) 
binom_prob <- do.call(rbind, lapply(1:n, function(i) dbinom(x = count, size = total, prob = theta[i], log = T)))
sum(matrixStats::colLogSumExps(binom_prob + beta_dens) - matrixStats::logSumExp(beta_dens))
