marginal_likelihood <- function(x, N, a, b) {
  # Binomial coefficient
  binom_coeff <- choose(N, x)
  
  # Beta function for posterior and prior
  beta_posterior <- beta(x + a, N - x + b)
  beta_prior <- beta(a, b)
  
  # Compute marginal likelihood
  marginal <- binom_coeff * beta_posterior / beta_prior
  return(marginal)
}

numerical_marginal_likelihood <- function(x, N, a, b, bins = 1000) {
  theta <- seq(0, 1, length.out = bins + 1) # Break interval [0, 1] into bins
  delta_theta <- 1 / bins # Bin width
  
  # Compute prior * likelihood for each theta
  integrand <- function(theta) {
    dbeta(theta, a, b) * dbinom(x, N, theta)
  }
  
  # Evaluate integrand at midpoints of bins
  midpoints <- (theta[-1] + theta[-length(theta)]) / 2
  integral <- sum(integrand(midpoints)) * delta_theta
  
  return(integral)
}

#####

x <- 173
N <- 1999
a <- 100
b <- 0.01

marginal_likelihood(x, N, a, b)
numerical_marginal_likelihood(x, N, a, b)
extraDistr::dbbinom(x, N, a, b)
