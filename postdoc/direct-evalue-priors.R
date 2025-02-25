# --------------------------------------------------
# 1) Simulate data with compound symmetry correlation
# --------------------------------------------------
library(MASS)   # for mvrnorm
library(cmdstanr)

p <- 8
n <- 10
rho <- 0.5

# Compound symmetry correlation matrix
corr_true <- matrix(rho, nrow=p, ncol=p)
diag(corr_true) <- 1

# Simulate data
X <- mvrnorm(n, mu = rep(0, p), Sigma = corr_true)

# ------------------------------------------
# 2) Stan model code with custom eigenvalue prior
# ------------------------------------------
stan_code <- '
data {
  int<lower=1> n;
  int<lower=1> p;
  matrix[n, p] x;
  int<lower=0, upper=1> incl_likelihood;
}

parameters {
  // Cholesky factor of a correlation matrix
  cholesky_factor_corr[p] L;
}

transformed parameters {
  // Build full correlation matrix
  matrix[p,p] R = multiply_lower_tri_self_transpose(L);
  // Compute its eigenvalues
  vector[p] lambda = eigenvalues_sym(R);
}

model {
  // *** 1) PRIOR on EIGENVALUES ***
  // For demonstration, let\'s say each eigenvalue ~ Normal(1, 0.5).
  // This is just an example. In practice, you can do any function of lambda.
  for(i in 1:p) {
    target += normal_lpdf(log(lambda[i]) | 1.0, 2.0);
  }

// *** 2) LIKELIHOOD (optional) ***
  if(incl_likelihood == 1) {
    for(i in 1:n) {
      target += multi_normal_lpdf(x[i] | rep_vector(0.0, p), R);
    }
  }
}

generated quantities {
  // Track the first off-diagonal correlation to see if it matches the data
  real offdiag = R[1,2];
}
'

# Compile the Stan model
mod <- stan_model(model_code = stan_code)

# ------------------------------------------
# 3) Fit model WITH likelihood
# ------------------------------------------
dat_with_lik <- list(n=n, p=p, x=X, incl_likelihood=1)
fit_with_lik  <- sampling(
  mod, data=dat_with_lik, 
  chains=2, iter=2000, warmup=500, 
  seed=42
)
print(fit_with_lik, pars=c("offdiag"), probs=c(0.1,0.5,0.9))

# ------------------------------------------
# 4) Fit model WITHOUT likelihood (PRIOR only)
# ------------------------------------------
dat_prior_only <- list(n=n, p=p, x=X, incl_likelihood=0)
fit_prior_only <- sampling(
  mod, data=dat_prior_only, 
  chains=2, iter=2000, warmup=500,
  seed=43
)
print(fit_prior_only, pars=c("offdiag"), probs=c(0.1,0.5,0.9))
